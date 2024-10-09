import math
from typing import Generator, Union
import pandas as pd
from pandas import Series, DataFrame

from logger import set_logger

logger = set_logger(name=__file__)

pd.options.mode.copy_on_write = True

CLONOTYPE_COLUMNS = ['v_call', 'j_call', 'junction']
V_ALIGN_COLUMN = 'v_sequence_alignment'
C_CALL_COLUMN = 'c_call'
J_CALL_COLUMN = 'j_call'
JUNCTION_COLUMN = 'junction'
COUNT_COLUMN = 'duplicate_count'
BASES_GAP = ['A', 'T', 'G', 'C', '']


class ClonotypeCounter:
    def __init__(self, v_call: str, j_call: str, junction: str, count: int, error_rate: float):
        self.junction = junction
        self.v_call = v_call
        self.j_call = j_call
        self.count = count
        self.parent = junction

        # probability to get single error according to Binomial distribution
        factor = min(0.1, error_rate * len(junction))
        self.factor = factor * (1.0 - factor)

        # 95% CI ~ confidence / sqrt(N), confidence intervals for proportion
        self.confidence = 1.96 * math.sqrt(self.factor * (1.0 - self.factor))

    def reassign_parent(
            self, new_v_call: str, new_j_call: str, new_junction: str, new_count: int, match_vj: bool = False
    ) -> None:
        if not match_vj or (new_v_call == self.v_call and new_j_call == self.j_call):
            if self.parent == self.junction and \
                    self.count / new_count < self.factor + self.confidence / math.sqrt(new_count):
                self.parent = new_junction

    def __repr__(self):
        return str(self.__dict__)


class ClonotypeCorrector:
    def __init__(self, top_c_call: bool, top_v_alignment_call: bool, error_rate: float = 0.001):
        self.error_rate = error_rate  # 0.001 - Phred Q30
        self.top_c_call = top_c_call
        self.top_v_alignment_call = top_v_alignment_call

    def correct_full(self, annotation: pd.DataFrame) -> pd.DataFrame:
        """
        Aggregates clonotypes and summarizes their duplicate count. Corrects errors in junction.
        Selects optimal C call (isotype) and top V alignment by coverage if requested.
        """
        annotation = annotation.reset_index(drop=True)
        aggregated_annotation = self.aggregate_clonotypes(annotation, CLONOTYPE_COLUMNS)
        fetched_annotation = self.fetch_clonotypes(aggregated_annotation)
        corrected_annotation = (self.correct_clonotypes(fetched_annotation)
                                .drop(columns=['count'])
                                .drop_duplicates())
        merged_annotation = aggregated_annotation.merge(corrected_annotation,
                                                        on=CLONOTYPE_COLUMNS,
                                                        how='left')
        full_corrected_annotation = (self.aggregate_clonotypes(merged_annotation, ['parent'])
                                     .drop(columns=['parent', 'factor']))
        return full_corrected_annotation

    def _aggregate_by_top_c_call(self, group: pd.DataFrame, clonotype: pd.Series) -> Series:
        c_calls = group.dropna(subset=C_CALL_COLUMN)
        if not c_calls.empty:
            c_calls = c_calls.groupby(C_CALL_COLUMN)[COUNT_COLUMN].sum()
            clonotype[C_CALL_COLUMN] = c_calls.nlargest(1).index[0]
        return clonotype

    def _aggregate_by_top_v_alignment_call(self, group: pd.DataFrame) -> DataFrame:
        v_alns = group.dropna(subset=V_ALIGN_COLUMN)
        if not v_alns.empty:
            v_alns['bases_count'] = v_alns[V_ALIGN_COLUMN].str.len() * v_alns[COUNT_COLUMN]
            return v_alns.nlargest(1, 'bases_count').drop('bases_count', axis=1)
        return group.nlargest(1, COUNT_COLUMN)

    def _aggregate_clonotypes_group(self, group: pd.DataFrame) -> DataFrame:
        """
        Selects annotation with the highest read count among all annotations for a given clonotype.
        Selects optimal C call (isotype) and top V alignment by coverage if requested.
        Optimal C call and V alignment are added as separate columns not included in AIRR format prefixed with 'best_'
        """
        clonotype = group.nlargest(1, COUNT_COLUMN)
        if self.top_v_alignment_call:
            v_alns = group.dropna(subset=V_ALIGN_COLUMN)
            if not v_alns.empty:
                v_alns['bases_count'] = v_alns[V_ALIGN_COLUMN].str.len() * v_alns[COUNT_COLUMN]
                best_v_aln = v_alns.nlargest(1, 'bases_count').drop('bases_count', axis=1).iloc[0]
                clonotype['best_v_call'] = best_v_aln['v_call']
                clonotype['best_v_sequence_alignment'] = best_v_aln['v_sequence_alignment']
                clonotype['best_v_sequence_alignment_aa'] = best_v_aln['v_sequence_alignment_aa']
                clonotype['best_v_germline_alignment'] = best_v_aln['v_germline_alignment']
                clonotype['best_v_germline_alignment_aa'] = best_v_aln['v_germline_alignment_aa']

        if self.top_c_call:
            c_calls = group.dropna(subset=C_CALL_COLUMN)
            if not c_calls.empty:
                c_calls = c_calls.groupby(C_CALL_COLUMN)[COUNT_COLUMN].sum()
                clonotype['best_c_call'] = c_calls.nlargest(1).index[0]

        return clonotype

    def aggregate_clonotypes(self, annotation: pd.DataFrame, grouping_columns: list[str]) -> Union[Series, DataFrame]:
        """Aggregates clonotypes by chosing optimal annotation and summarizing duplicate count
        across all annotations for a given clonotype."""
        annotation = annotation.reset_index(drop=True)
        clonotype_groups = annotation.groupby(grouping_columns)

        duplicate_count = clonotype_groups[COUNT_COLUMN].transform('sum')
        aggregated_groups = [self._aggregate_clonotypes_group(group) for _, group in clonotype_groups]

        if not aggregated_groups:
            return annotation

        aggregated_annotation = pd.concat(aggregated_groups)
        aggregated_annotation[COUNT_COLUMN] = duplicate_count
        aggregated_annotation.sort_values(by=COUNT_COLUMN, ascending=False, inplace=True)

        logger.info(f'Filtered out {len(annotation) - len(aggregated_annotation)} clones while aggregation.')

        return aggregated_annotation

    def fetch_clonotypes(self, annotation: pd.DataFrame) -> pd.DataFrame:
        fetched_annotation = (
            annotation[annotation[JUNCTION_COLUMN].values != '']
            .groupby(CLONOTYPE_COLUMNS)[COUNT_COLUMN]
            .sum()
            .reset_index()
            .rename(columns={'sum': COUNT_COLUMN})
            .sort_values(by=COUNT_COLUMN, ascending=False)
        )
        logger.info(f'Filtered out {annotation.shape[0] - fetched_annotation.shape[0]} clones while fetching.')
        return fetched_annotation

    def correct_clonotypes(self, annotation: pd.DataFrame) -> pd.DataFrame:
        """
        Runs frequency-based error correction, assigns all clonotypes
        that differ by a mismatch in junction with a parent clonotype.

        Will assign a clonotype to another one
        if the smallest (child) has 'threshold' times less duplicate count compared to largest (parent).
        """
        clonotype_count_list, clonotype_count_dict = self._make_counters(annotation)
        corrected_annotation = self._update_counters_inplace(clonotype_count_list=clonotype_count_list,
                                                             clonotype_count_dict=clonotype_count_dict)
        logger.info(f'Filtered out {annotation.shape[0] - corrected_annotation.shape[0]} clones while correcting.')
        return corrected_annotation

    def _make_counters(self, annotation: pd.DataFrame) -> tuple[list, dict]:
        """Aux. routine for initializing clonotype duplicate counters for frequency-based error correction"""
        annotation = annotation.sort_values(by=COUNT_COLUMN, ascending=False)

        clonotype_count_dict = {}
        clonotype_count_list = [] # note that list is sorted from largest to smallest, needed for _update_counters_inplace

        for clonotype, count in zip(annotation[CLONOTYPE_COLUMNS].values, annotation[COUNT_COLUMN].values):
            counter = ClonotypeCounter(*clonotype, count, error_rate=self.error_rate)
            clonotype_count_dict[counter.junction] = clonotype_count_dict.get(counter.junction, []) + [counter]
            clonotype_count_list.append(counter)

        return clonotype_count_list, clonotype_count_dict

    def _update_counters_inplace(self, clonotype_count_list: list, clonotype_count_dict: dict) -> pd.DataFrame:
        for counter in clonotype_count_list:
            for junction_variant in self._get_variants(counter.junction): # we search using current junction
                for counter_to_update in clonotype_count_dict.get(junction_variant, []): # we assign parent's junction as parent
                    counter_to_update.reassign_parent(counter.v_call, counter.j_call, counter.parent, counter.count)
        return (pd.DataFrame
                .from_records([count.__dict__ for count in clonotype_count_list])
                .rename(columns={'seq': JUNCTION_COLUMN}))

    def _get_variants(self, sequence: str) -> Generator:
        """Generates variants of the sequence"""
        for i, bp in enumerate(sequence):
            yield from self._replace_base(sequence, i, bp)
            yield from self._insert_base(sequence, i)

    def _replace_base(self, sequence: str, index: int, base: str) -> Generator:
        """Yields sequences with the base at the given index replaced"""
        for bp_new in BASES_GAP:
            if base != bp_new:
                yield sequence[:index] + bp_new + sequence[(index + 1):]

    def _insert_base(self, sequence: str, index: int) -> Generator:
        """Yields sequences with new bases inserted after the given index"""
        for bp_new in BASES_GAP:
            if bp_new:
                yield sequence[:index + 1] + bp_new + sequence[(index + 1):]
