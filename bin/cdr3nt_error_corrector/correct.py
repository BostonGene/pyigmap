import pandas as pd
import math
from typing import Generator
from logger import set_logger

logger = set_logger(name=__file__)

pd.options.mode.copy_on_write = True


class ClonotypeCounter:
    def __init__(self, v_call: str, j_call: str, junction: str, count: int, error_rate: float):
        self.junction = junction
        self.v_call = v_call
        self.j_call = j_call
        self.count = count
        self.parent = junction
        # probabiliy to get single error according to Binomial distr
        self.factor = min(0.1, error_rate * len(junction) )
        self.factor = self.factor * (1.0 - self.factor)
        # 95% CI ~ confidence / sqrt(N), confidence intervals for proportion
        self.confidence = 1.96 * math.sqrt(self.factor * (1.0 - self.factor)) 

    def reassign_parent(self, v_call: str, j_call: str, junction: str, count: int, matchVJ: bool = False):
        if not matchVJ or (v_call == self.v_call and j_call == self.j_call):
            if self.parent == self.junction and self.count / count < self.factor + self.confidence / math.sqrt(count):
                self.parent = junction

    def __repr__(self):
        return str(self.__dict__)


class ClonotypeCorrector:
    CLONOTYPE_COLUMNS = ['v_call', 'j_call', 'junction']
    V_ALIGN_COLUMN = 'v_sequence_alignment'
    C_CALL_COLUMN = 'c_call'
    J_CALL_COLUMN = 'j_call'
    JUNCTION_COLUMN = 'junction'
    COUNT_COLUMN = 'duplicate_count'
    BASES_GAP = ['A', 'T', 'G', 'C', '']

    def __init__(self, top_c_call: bool, top_v_alignment_call: bool, error_rate=0.001):
        self.error_rate = error_rate # 0.001 - Phred Q30
        self.top_c_call = top_c_call
        self.top_v_alignment_call = top_v_alignment_call

    def _correct_c_call(self, c_call: str, j_call: str) -> str:
        if j_call.startswith("TR"):
            c_call = j_call[0:3] + "C"
            if c_call == "TRBC":
                c_call = c_call + j_call[4]
        return c_call

    def correct_full(self, annotation: pd.DataFrame) -> pd.DataFrame:
        """
        Aggregates clonotypes and summarizes their duplicate count. Corrects errors in junction. 
        Selects optimal C call (isotype) and top V alignment by coverage if requested.
        """
        annotation = annotation.reset_index(drop=True)
        annotation[self.C_CALL_COLUMN] = annotation.apply(lambda x: self._correct_c_call(x[self.C_CALL_COLUMN], x[self.J_CALL_COLUMN]), 
                                                          axis=1) # simple fix for lack of TRxC genes
        aggregated_annotation = self.aggregate_clonotypes(annotation, self.CLONOTYPE_COLUMNS)
        fetched_annotation = self.fetch_clonotypes(aggregated_annotation)
        corrected_annotation = (self.correct_clonotypes(fetched_annotation)
                                .drop(columns=['count'])
                                .drop_duplicates())
        merged_annotation = aggregated_annotation.merge(corrected_annotation,
                                                        on=self.CLONOTYPE_COLUMNS,
                                                        how='left')
        full_corrected_annotation = (self.aggregate_clonotypes(merged_annotation, ['parent'])
                                     .drop(columns=['parent', 'factor']))
        return full_corrected_annotation

    def _aggregate_clonotypes_group(self, group: pd.DataFrame) -> pd.Series:
        """
        Selects annotation with highest read count among all annotations for a given clonotype.
        Selects optimal C call (isotype) and top V alignment by coverage if requested.         
        """
        if self.top_v_alignment_call:
            v_alns = group.dropna(subset=self.V_ALIGN_COLUMN)            
            if not v_alns.empty:
                v_alns['bases_count'] = v_alns[self.V_ALIGN_COLUMN].str.len() * v_alns[self.COUNT_COLUMN]
                clonotype = v_alns.nlargest(1, 'bases_count').drop('bases_count', axis=1)
            else:
                clonotype = group.nlargest(1, self.COUNT_COLUMN)
        else:
            clonotype = group.nlargest(1, self.COUNT_COLUMN)

        if self.top_c_call:
            c_calls = group.dropna(subset=self.C_CALL_COLUMN)
            if not c_calls.empty:
                c_calls = c_calls.groupby(self.C_CALL_COLUMN)[self.COUNT_COLUMN].sum()
                clonotype[self.C_CALL_COLUMN] = c_calls.nlargest(1).index[0]
                # TODO: also select other C-alignment-related columns

        return clonotype

    def aggregate_clonotypes(self, annotation: pd.DataFrame, grouping_columns: list[str]) -> pd.DataFrame:
        """
        Aggregates clonotypes by chosing optimal annotation and summarizing duplicate count across all annotations for a given clonotype.
        """
        annotation = annotation.reset_index(drop=True)
        clonotype_groups = annotation.groupby(grouping_columns)
        duplicate_count = clonotype_groups[self.COUNT_COLUMN].transform('sum')
        aggregated_groups = [self._aggregate_clonotypes_group(group) for _, group in clonotype_groups]
        if not aggregated_groups:
            return annotation
        aggregated_annotation = pd.concat(aggregated_groups)
        aggregated_annotation[self.COUNT_COLUMN] = duplicate_count
        aggregated_annotation.sort_values(by=self.COUNT_COLUMN, ascending=False, inplace=True)
        logger.info(f'Filtered out {len(annotation) - len(aggregated_annotation)} clones while aggregation.')
        return aggregated_annotation

    def fetch_clonotypes(self, annotation: pd.DataFrame) -> pd.DataFrame:
        """
        Gets clontypes defined as a tuple with junction, V and J calls
        """
        fetched_annotation = (annotation[annotation[self.JUNCTION_COLUMN].values != '']
                              .groupby(self.CLONOTYPE_COLUMNS)[self.COUNT_COLUMN]
                              .sum()
                              .reset_index()
                              .rename(columns={'sum': self.COUNT_COLUMN})
                              .sort_values(by=self.COUNT_COLUMN, ascending=False))
        logger.info(f'Filtered out {annotation.shape[0] - fetched_annotation.shape[0]} clones while fetching.')
        return fetched_annotation

    def correct_clonotypes(self, annotation: pd.DataFrame):
        """
        Runs frequency-based error correction, assigns all clonotypes that differ by a mismatch in junction with a parent clonotype.
        Will assign a clonotype to another one if the smallest (child) has 'threshold' times less duplicate count compared to largest (parent).
        """
        clonotype_count_list, clonotype_count_dict = self._make_counters(annotation)
        corrected_annotation = self._update_counters_inplace(clonotype_count_list=clonotype_count_list,
                                                             clonotype_count_dict=clonotype_count_dict)
        logger.info(f'Filtered out {annotation.shape[0] - corrected_annotation.shape[0]} clones while correcting.')
        return corrected_annotation

    def _make_counters(self, annotation: pd.DataFrame):
        """
        Aux. routine for initializing clonotype duplicate counters for frequency-based error correction
        """
        annotation = annotation.sort_values(by=self.COUNT_COLUMN, ascending=False)

        clonotype_count_dict = {}
        clonotype_count_list = []

        for clonotype, count in zip(annotation[self.CLONOTYPE_COLUMNS].values, annotation[self.COUNT_COLUMN].values):
            counter = ClonotypeCounter(*clonotype, count, error_rate=self.error_rate)
            clonotype_count_dict[counter.junction] = clonotype_count_dict.get(counter.junction, []) + [counter]
            clonotype_count_list.append(counter)

        return clonotype_count_list, clonotype_count_dict

    def _update_counters_inplace(self, clonotype_count_list: list, clonotype_count_dict: dict) -> pd.DataFrame:
        """
        Try to update counters and reassign parent for all clonotype variants that differ by a single nucleotide mismatch in junction.
        """
        for counter in clonotype_count_list:
            for junction_variant in self._get_variants(counter.junction):
                for counter_to_update in clonotype_count_dict.get(junction_variant, []):
                    counter_to_update.reassign_parent(counter.v_call, counter.j_call, counter.junction, counter.count)
        return (pd.DataFrame
                .from_records([count.__dict__ for count in clonotype_count_list])
                .rename(columns={'seq': self.JUNCTION_COLUMN}))

    def _get_variants(self, sequence: str) -> Generator:
        """
        Generates variants of the sequence.

        :param sequence: A sequence string.

        return: Sequence with one replaced nucleotide.
        """
        for i, bp in enumerate(sequence):
            for bp_new in self.BASES_GAP:
                if bp != bp_new:
                    yield sequence[:i] + bp_new + sequence[(i + 1):]
                if bp_new:
                    yield sequence[:i + 1] + bp_new + sequence[(i + 1):]
