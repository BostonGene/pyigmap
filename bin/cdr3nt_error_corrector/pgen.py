from multiprocessing import Pool
import os
from typing import Union

import olga.load_model as load_model
import olga.generation_probability as pgen
from olga.generation_probability import GenerationProbabilityVDJ, GenerationProbabilityVJ

from logger import set_logger

logger = set_logger(name=__file__)

LOCUS_GLOSSARY = {
    'TRA': 'T_alpha',
    'TRB': 'T_beta',
    'IGH': 'B_heavy',
    'IGK': 'B_kappa',
    'IGL': 'B_lambda'
}

VDJ_LOCI = {"IGH", "TRB"}  # These loci include Variable (V), Diversity (D), and Joining (J) gene segments


class PgenModel:
    """Wrapper of the OLGA tool https://github.com/statbiophys/OLGA"""

    def __init__(self, olga_models_dir: str, locus: str, specie: str = 'human'):
        self.locus = locus[:3]
        self.specie = specie
        self.olga_models_dir = olga_models_dir
        self.model = self.get_olga_model()

    def get_pgen(self, cdr3_aa: list[str]) -> list[float]:
        """Calculates the generation probability (Pgen) of the sequence in parallel"""
        logger.info(f'Detecting spurious rearrangements in {self.locus} via OLGA tool...')

        with Pool(processes=os.cpu_count()) as pool:
            pgen_values = pool.map(self.calculate_pgen, cdr3_aa)

        logger.info(f'Spurious rearrangements detection in {self.locus} has been done.')

        return pgen_values

    def calculate_pgen(self, cdr3_aa: str) -> Union[float, None]:
        """Calculates the generation probability (Pgen) of the sequence"""
        if self.model:
            probability = self.model.compute_aa_CDR3_pgen(cdr3_aa)
            return probability
        # None for CDR3 with locus that does not exist in OLGA model (for example: TRG, TRD)
        # https://github.com/statbiophys/OLGA/tree/master/olga/default_models
        return None

    def _get_data_class_for_locus(self) -> Union[load_model.GenomicDataVDJ, load_model.GenomicDataVJ]:
        return load_model.GenomicDataVDJ() if self.locus in VDJ_LOCI else load_model.GenomicDataVJ()

    def _get_model_for_locus(self) -> Union[load_model.GenerativeModelVDJ, load_model.GenerativeModelVJ]:
        return load_model.GenerativeModelVDJ() if self.locus in VDJ_LOCI \
            else load_model.GenerativeModelVJ()

    def _get_pgen_for_locus(self) -> Union[GenerationProbabilityVDJ, GenerationProbabilityVJ]:
        return pgen.GenerationProbabilityVDJ if self.locus in VDJ_LOCI else pgen.GenerationProbabilityVJ

    def _get_igor_model_files(self, locus_folder_name: str) -> tuple[str, str, str, str]:
        model_params_file = os.path.join(self.olga_models_dir, locus_folder_name, 'model_params.txt')
        model_marginals_file = os.path.join(self.olga_models_dir, locus_folder_name, 'model_marginals.txt')
        v_anchor_pos_file = os.path.join(self.olga_models_dir, locus_folder_name, 'V_gene_CDR3_anchors.csv')
        j_anchor_pos_file = os.path.join(self.olga_models_dir, locus_folder_name, 'J_gene_CDR3_anchors.csv')

        return model_params_file, model_marginals_file, v_anchor_pos_file, j_anchor_pos_file

    def get_olga_model(self) -> Union[GenerationProbabilityVDJ, GenerationProbabilityVJ, None]:
        if LOCUS_GLOSSARY.get(self.locus):
            locus_glossary = LOCUS_GLOSSARY.get(self.locus)
            locus_folder_name = f'{self.specie}_{locus_glossary}'

            model_params_file, model_marginals_file, v_anchor_pos_file, j_anchor_pos_file \
                = self._get_igor_model_files(locus_folder_name)

            generative_model = self._get_model_for_locus()
            generative_model.load_and_process_igor_model(model_marginals_file)

            genomic_data = self._get_data_class_for_locus()
            genomic_data.load_igor_genomic_data(model_params_file, v_anchor_pos_file, j_anchor_pos_file)

            probability_generator = self._get_pgen_for_locus()

            return probability_generator(generative_model, genomic_data)

        logger.warning(f"{self.locus} does not exist in OLGA models, pgen calculation in this locus will be skipped...")
        # None for locus that does not exist in OLGA model (for example: TRG, TRD)
        # https://github.com/statbiophys/OLGA/tree/master/olga/default_models
        return None
