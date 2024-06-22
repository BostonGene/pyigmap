import pandas as pd
from typing import Generator
from logger import set_logger

logger = set_logger(name=__file__)

pd.options.mode.copy_on_write = True


class ClonotypeCounter:
    def __init__(self, v_call: str, j_call: str, junction: str, count: int, factor: float):
        self.junction = junction
        self.v_call = v_call
        self.j_call = j_call
        self.count = count
        self.parent = junction
        self.factor = factor

    def reassign_parent(self, v_call: str, j_call: str, junction: str, count: int, matchVJ: bool = False):
        if not matchVJ or (v_call == self.v_call and j_call == self.j_call):
            if self.parent == self.junction and (self.count + 1) / (count + 1) < self.factor:
                self.parent = junction

    def __repr__(self):
        return str(self.__dict__)


class ClonotypeCorrector:
    CLONOTYPE_COLUMNS = ['v_call', 'j_call', 'junction']
    V_ALIGN_COLUMNS = ['v_sequence_alignment', 'v_germline_alignment', 'v_cigar']
    C_CALL_COLUMN = 'c_call'
    JUNCTION_COLUMN = 'junction'
    COUNT_COLUMN = 'duplicate_count'
    BASES_GAP = ['A', 'T', 'G', 'C', '']

    def __init__(self, top_c_call: bool, top_v_alignment_call: bool, collapse_factor=None):
        self.factor = collapse_factor or 0.05
        self.top_c_call = top_c_call
        self.top_v_alignment_call = top_v_alignment_call

    def correct_full(self, annotation: pd.DataFrame) -> pd.DataFrame:
        annotation = annotation.reset_index(drop=True)

        aggregated_annotation = self.aggregate_full(annotation, self.CLONOTYPE_COLUMNS)

        fetched_annotation = self.fetch_clonotypes(aggregated_annotation)

        corrected_annotation = (self.correct_clonotypes(fetched_annotation)
                                .drop(columns=['count'])
                                .drop_duplicates())

        merged_annotation = aggregated_annotation.merge(corrected_annotation,
                                                        on=self.CLONOTYPE_COLUMNS,
                                                        how='left')

        full_corrected_annotation = (self.aggregate_full(merged_annotation, ['parent'])
                                     .drop(columns=['parent', 'factor']))

        return full_corrected_annotation

    def aggregate_full(self, annotation: pd.DataFrame, grouping_columns: list[str]) -> pd.DataFrame:
        """Runs full clonotypes aggregation using c_call and v alignment columns"""
        aggregated_annotation_c_call = self.aggregate_clonotypes(annotation[~annotation[self.C_CALL_COLUMN].isna()],
                                                                 grouping_columns + [self.C_CALL_COLUMN])
        aggregated_annotation_no_c_call = self.aggregate_clonotypes(annotation[annotation[self.C_CALL_COLUMN].isna()],
                                                                    grouping_columns)

        aggregated_annotation = pd.concat([aggregated_annotation_no_c_call, aggregated_annotation_c_call])

        return aggregated_annotation

    def _most_weighted(self, group: pd.DataFrame) -> pd.DataFrame:
        """Returns clonotypes with the biggest weight (count)"""
        group['count'] = group[self.COUNT_COLUMN].sum()
        group['max_count'] = group.groupby(self.CLONOTYPE_COLUMNS)['count'].transform('max')
        return group[group['count'] == group['max_count']].drop(columns=['count', 'max_count'])

    def _most_frequent(self, group: pd.DataFrame, columns: list[str]):
        """Returns clonotype with the most frequent values in selected columns"""
        group_filtered = group.dropna(subset=columns)

        if group_filtered.empty:
            return group

        value_counts = group_filtered[columns].value_counts().reset_index(name='count')
        most_frequent_row = value_counts.loc[value_counts['count'].idxmax()]
        most_frequent_values = most_frequent_row[columns].to_dict()

        clonotypes_with_most_frequent_values = group_filtered[
            (group_filtered[columns] == pd.Series(most_frequent_values)).all(axis=1)
        ]

        return clonotypes_with_most_frequent_values

    def _aggregate_clonotypes_group(self, group: pd.DataFrame) -> pd.Series:
        """Select the most weighted clonotype with the most frequent values"""
        weighted_group = self._most_weighted(group)

        top_clonotypes = weighted_group

        if self.top_c_call:
            top_clonotypes = self._most_frequent(weighted_group, [self.C_CALL_COLUMN])

        if self.top_v_alignment_call:
            top_clonotypes = self._most_frequent(top_clonotypes, self.V_ALIGN_COLUMNS)

        return top_clonotypes.head(1)

    def aggregate_clonotypes(self, annotation: pd.DataFrame, grouping_columns: list[str]) -> pd.DataFrame:
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
        fetched_annotation = (annotation[annotation[self.JUNCTION_COLUMN].values != '']
                              .groupby(self.CLONOTYPE_COLUMNS)[self.COUNT_COLUMN]
                              .sum()
                              .reset_index()
                              .rename(columns={'sum': self.COUNT_COLUMN})
                              .sort_values(by=self.COUNT_COLUMN, ascending=False))
        logger.info(f'Filtered out {annotation.shape[0] - fetched_annotation.shape[0]} clones while fetching.')
        return fetched_annotation

    def correct_clonotypes(self, annotation: pd.DataFrame):
        clonotype_count_list, clonotype_count_dict = self._make_counters(annotation)
        corrected_annotation = self._update_counters_inplace(clonotype_count_list=clonotype_count_list,
                                                             clonotype_count_dict=clonotype_count_dict)
        logger.info(f'Filtered out {annotation.shape[0] - corrected_annotation.shape[0]} clones while correcting.')
        return corrected_annotation

    def _make_counters(self, annotation: pd.DataFrame):
        annotation = annotation.sort_values(by=self.COUNT_COLUMN, ascending=False)

        clonotype_count_dict = {}
        clonotype_count_list = []

        for clonotype, count in zip(annotation[self.CLONOTYPE_COLUMNS].values, annotation[self.COUNT_COLUMN].values):
            counter = ClonotypeCounter(*clonotype, count, factor=self.factor)
            clonotype_count_dict[counter.junction] = clonotype_count_dict.get(counter.junction, []) + [counter]
            clonotype_count_list.append(counter)

        return clonotype_count_list, clonotype_count_dict

    def _update_counters_inplace(self, clonotype_count_list: list, clonotype_count_dict: dict) -> pd.DataFrame:
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
