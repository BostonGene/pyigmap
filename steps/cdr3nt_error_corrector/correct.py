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
    JUNCTION_COLUMN = 'junction'
    COUNT_COLUMN = 'duplicate_count'
    BASES_GAP = ['A', 'T', 'G', 'C', '']

    def __init__(self, use_v_align: bool, collapse_factor=None):
        self.factor = collapse_factor or 0.05
        self.use_v_align = use_v_align

    def correct_full(self, annotation: pd.DataFrame) -> pd.DataFrame:
        annotation = annotation.reset_index(drop=True)

        aggregated_annotation = self.aggregate_clonotypes(annotation, self.CLONOTYPE_COLUMNS)

        fetched_annotation = self.fetch_clonotypes(aggregated_annotation)

        corrected_annotation = (self.correct_clonotypes(fetched_annotation)
                                .drop(columns=['count'])
                                .drop_duplicates())

        merged_annotation = aggregated_annotation.merge(corrected_annotation,
                                                        on=self.CLONOTYPE_COLUMNS,
                                                        how='left')

        full_corrected_annotation = (self.aggregate_clonotypes(merged_annotation, grouping_columns=['parent'])
                                     .drop(columns=['parent', 'factor']))

        return full_corrected_annotation

    def _most_weighted(self, group: pd.DataFrame) -> pd.DataFrame:
        """Returns clonotypes with the biggest weight (count)"""
        group['count'] = group[self.COUNT_COLUMN].sum()
        group['max_count'] = group.groupby(self.CLONOTYPE_COLUMNS)['count'].transform('max')
        return group[group['count'] == group['max_count']]

    @staticmethod
    def _most_frequent_values(df, columns):
        most_frequent_values = {}

        df_without_na = df.dropna(subset=columns)

        if not len(df_without_na):
            return df

        for col in columns:
            most_freq_value = df_without_na[col].value_counts().idxmax()
            most_frequent_values[col] = most_freq_value

        condition = df_without_na[columns].eq(pd.Series(most_frequent_values))
        most_frequent_row_df = df_without_na[condition.all(axis=1)]

        return most_frequent_row_df.head(1)

    def _aggregate_clonotypes_group(self, group: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
        """Returns the most weighted clonotype and with the most frequent values in selected columns"""

        group_weighted = self._most_weighted(group)

        columns_to_drop = list(set(columns).union(set(self.CLONOTYPE_COLUMNS)))

        group_weighted.drop(
            columns=group_weighted.columns.difference(columns_to_drop),
            inplace=True, axis=1
        )

        if group_weighted[columns].isna().all().all():
            return group_weighted

        # most_frequent = group_weighted[columns].mode(dropna=True).iloc[0]
        #
        # mask = (group_weighted[columns] == most_frequent).all(axis=1)

        # return group_weighted[mask].head(1)

        return ClonotypeCorrector._most_frequent_values(group_weighted, columns)

    def aggregate_clonotypes(self, annotation: pd.DataFrame, grouping_columns: list) -> pd.DataFrame:
        annotation = annotation.reset_index(drop=True)

        clonotype_groups = annotation.groupby(grouping_columns)

        duplicate_count = clonotype_groups[self.COUNT_COLUMN].transform('sum')

        aggregated_c_call_groups = pd.concat([self._aggregate_clonotypes_group(group, ['c_call'])
                                              for _, group in clonotype_groups])
        columns_to_drop = ['c_call']

        if self.use_v_align:
            aggregated_v_align_groups = pd.concat([self._aggregate_clonotypes_group(group, self.V_ALIGN_COLUMNS)
                                                   for _, group in clonotype_groups])
            columns_to_drop.append(self.V_ALIGN_COLUMNS)

        annotation = annotation.drop(columns=columns_to_drop)
        clonotype_groups = annotation.groupby(grouping_columns)
        aggregated_annotation = pd.concat([self._aggregate_clonotypes_group(group, annotation.columns)
                                           for _, group in clonotype_groups])

        aggregated_annotation = aggregated_annotation.merge(aggregated_c_call_groups, how='left',
                                                            on=self.CLONOTYPE_COLUMNS)

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
