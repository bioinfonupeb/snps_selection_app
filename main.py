import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import List, Tuple
import logging
from datetime import datetime


class VariantLine:
    def __init__(self, line: str, samples_idx: List[int]) -> None:
        self.line = line.strip()
        if len(self.line) == 0:
            raise ValueError('The line must not be empty.')
        self.samples_idx = samples_idx
        if len(self.samples_idx) != 4:
            raise ValueError('The samples index must have 4 elements.')
        self.extract_info()

    def extract_info(self):
        fields = self.line.split('\t')
        samples_idx = self.samples_idx
        self.allele_ref = fields[3]
        self.allele_alt = fields[4]
        self.parental_sup_info = [x for x in fields[samples_idx[0]].split(':')]
        self.parental_inf_info = [x for x in fields[samples_idx[1]].split(':')]
        self.pool_sup_info = [x for x in fields[samples_idx[2]].split(':')]
        self.pool_rnd_info = [x for x in fields[samples_idx[3]].split(':')]
        self.sample_info = [self.parental_sup_info, self.parental_inf_info, self.pool_sup_info, self.pool_rnd_info]


class VCFProcessor:
    def __init__(self, gui_handler) -> None:
        self.gui_handler = gui_handler
        self.logger = self.gui_handler.logger
        self.vcf_file_path = None
        self.metadata_lines = []
        self.original_lines = []
        self.corrected_lines = []
        self.variant_lines = []
        self.header_fields = []
        self.header_idx = {}
        self.sample_columns = []  # Lista das colunas de amostras detectadas

        self.bug_lines = list()
        self.multiple_alleles_lines = list()
        self.not_single_genotype_lines = list()
        self.equal_genotype_lines = list()
        self.low_reads_lines = list()
        self.insufficient_diff_lines = list()

        # Initialize indices
        self.parental_sup_idx = None
        self.parental_inf_idx = None
        self.pool_sup_idx = None
        self.pool_rnd_idx = None
        self.alt_idx = None
        self.ref_idx = None

    def set_file_path(self, file_path: str) -> None:
        self.vcf_file_path = file_path
        self.logger.debug(f'vcf file path: {self.vcf_file_path}')

    def load_vcf_file(self) -> bool:
        '''Carrega e parseia o arquivo VCF'''
        try:
            with open(self.vcf_file_path, 'r') as f:
                lines = f.readlines()
                self.metadata_lines = [line for line in lines if line.startswith('#')]
                self.original_lines = [line for line in lines if not line.startswith('#')]
                self.variant_lines = self.original_lines.copy()

                # Encontrar a linha de header
                for line in self.metadata_lines:
                    if line.startswith('#CHROM'):
                        self.parse_header(line)
                        break

                self.gui_handler.log(f'Arquivo VCF carregado: {len(self.variant_lines)} variantes encontradas')
                # self.gui_handler.log(f'Colunas de amostras detectadas: {self.sample_columns}')
                return True
        except Exception as e:
            self.gui_handler.log(f'Erro ao ler arquivo VCF: {str(e)}')
            return False

    def parse_header(self, header_line: str):
        '''Parseia a linha de header para detectar colunas de amostras'''
        header_line = header_line.strip()
        if header_line.startswith('#'):
            header_line = header_line[1:]  # Remove o #

        self.header_fields = header_line.split('\t')
        self.header_idx = {field: idx for idx, field in enumerate(self.header_fields)}

        # As colunas de amostras são as últimas 4 colunas (após FORMAT)
        format_idx = self.header_idx.get('FORMAT')
        if format_idx is not None and len(self.header_fields) >= format_idx + 5:
            self.sample_columns = self.header_fields[format_idx + 1:format_idx + 5]
        else:
            # Fallback: pegar as últimas 4 colunas
            self.sample_columns = self.header_fields[-4:]

        # Configurar índices básicos
        self.alt_idx = self.header_idx.get('ALT')
        self.ref_idx = self.header_idx.get('REF')

    def set_sample_columns(self, parental_sup: str, parental_inf: str, pool_sup: str, pool_rnd: str) -> bool:
        '''Configura as colunas das amostras baseado na seleção do usuário'''
        try:
            self.parental_sup_idx = self.header_idx[parental_sup]
            self.parental_inf_idx = self.header_idx[parental_inf]
            self.pool_sup_idx = self.header_idx[pool_sup]
            self.pool_rnd_idx = self.header_idx[pool_rnd]

            # self.gui_handler.log(f'Colunas configuradas: P-Sup={parental_sup}({self.parental_sup_idx}), '
            #                  f'P-Inf={parental_inf}({self.parental_inf_idx}), '
            #                  f'Pool-Sup={pool_sup}({self.pool_sup_idx}), '
            #                  f'Pool-Rnd={pool_rnd}({self.pool_rnd_idx})')
            return True
        except Exception as e:
            self.gui_handler.log(f'Erro ao configurar colunas: {str(e)}')
            return False

    def select_most_freq_genotype(self, sample_info: str, nuc_ref: str, nuc_alt: str, verbose: bool = False) -> str:
        '''
        Method to select the most frequent genotype in the sample info field and correct it to 0 or 1
        :param sample_info:
        :param nuc_ref:
        :param nuc_alt:
        :param verbose:
        :return:
        '''
        fields = sample_info.split(':')
        valids = ['0', '1', '.']
        genotypes = fields[0].split('/')
        flag = True

        for genotype in genotypes:
            if genotype not in valids:
                if verbose:
                    self.logger.info(f'Invalid genotype: {genotype}')
                flag = False
                break

        if flag:
            nuc, counts = self.get_major_nuc(fields[4])
            if nuc == nuc_ref:
                fields[0] = '0'
            elif nuc == nuc_alt:
                fields[0] = '1'
            else:
                if verbose:
                    self.logger.info(f'INVALID LINE: {fields[0]}, {nuc}, {nuc_ref}, {nuc_alt}, {counts}')
                return False

        corrected_samples_info = ':'.join(fields)
        return corrected_samples_info

    def correct_genotypes(self, verbose: bool = True) -> list[str]:
        '''
        Method to correct the genotypes of the variant lines
        :param verbose:
        :return:
        '''
        selected_rows = list()
        for i, row in enumerate(self.variant_lines):
            flag = True
            fields = row.split('\t')
            # print(f'{i} Fields: {fields}')
            nuc_ref = fields[self.ref_idx]
            nuc_alt = fields[self.alt_idx]
            # print(f'{i} Nuc Ref: {nuc_ref}, Nuc Alt: {nuc_alt}')
            if len(nuc_ref) == 1 and len(nuc_alt) == 1:
                for j in range(-4, -2):
                    sample = fields[j]
                    corrected = self.select_most_freq_genotype(sample, nuc_ref, nuc_alt, verbose=verbose)
                    # print(f'{i} Sample: {sample}, Corrected: {corrected}')
                    if corrected:
                        fields[j] = corrected
                    else:
                        flag = False
            else:
                if len(nuc_ref.split(',')) > 1 or len(nuc_alt.split(',')) > 1:
                    # self.logger.info(f'Variant with multiple alleles: {nuc_ref} - {nuc_alt}')
                    self.multiple_alleles_lines.append(row)
                else:
                    self.not_single_genotype_lines.append(row)
                continue

            new_line = '\t'.join([str(x) for x in fields])

            parental_inf_info = fields[self.parental_inf_idx].split(':')
            parental_sup_info = fields[self.parental_sup_idx].split(':')
            pool_sup_info = fields[self.pool_sup_idx].split(':')
            pool_rnd_info = fields[self.pool_rnd_idx].split(':')

            # Verify if parental_sup genotype and parental_inf genotype are equal
            genotype_sup = parental_sup_info[0]
            genotype_inf = parental_inf_info[0]

            if genotype_sup == genotype_inf:
                self.equal_genotype_lines.append(new_line)
                continue

            # Verify if allele counts are less than 3
            samples_fields = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]

            min_flag = False
            for sample in samples_fields:
                counts = self.get_counts(sample)
                if sum(counts) <= 3:
                    min_flag = True

            if min_flag:
                self.low_reads_lines.append(new_line)
                continue

            # Verify if the difference in the number of reads is less than 1
            parental_sup_counts = self.get_counts(parental_sup_info)
            parental_inf_counts = self.get_counts(parental_inf_info)
            parental_sup_ordered = sorted(parental_sup_counts, reverse=True)
            parental_inf_ordered = sorted(parental_inf_counts, reverse=True)

            min_diff = 0.2
            diff_perc_sup = abs(parental_sup_ordered[0] - parental_sup_ordered[1]) / sum(parental_sup_counts)
            diff_perc_inf = abs(parental_inf_ordered[0] - parental_inf_ordered[1]) / sum(parental_inf_counts)

            if diff_perc_sup < min_diff or diff_perc_inf < min_diff:
                self.insufficient_diff_lines.append(new_line)
                continue

            if flag:
                selected_rows.append(new_line)
            else:
                self.bug_lines.append(new_line)

        self.corrected_lines = selected_rows
        self.variant_lines = self.corrected_lines.copy()

        if verbose:
            self.gui_handler.log(f"Corrected genotypes successfully")
            self.gui_handler.log(f'Number of original lines: {len(self.original_lines)}')
            self.gui_handler.log(f'Number of corrected lines: {len(self.variant_lines)}')
            self.gui_handler.log(f'Number of not single genotype lines: {len(self.not_single_genotype_lines)}')
            self.gui_handler.log(f'Number of multiple alleles lines: {len(self.multiple_alleles_lines)}')
            self.gui_handler.log(f'Number of equal genotype lines: {len(self.equal_genotype_lines)}')
            self.gui_handler.log(f'Number of low reads lines: {len(self.low_reads_lines)}')
            self.gui_handler.log(f'Number of insufficient diff lines: {len(self.insufficient_diff_lines)}')
            self.gui_handler.log(f'Number of bug lines: {len(self.bug_lines)}')

        return selected_rows

    def get_info_fields(self, variant_line: str) -> Tuple[str, str, List[List[str]]]:
        '''Extrai informações dos campos da variante'''
        fields = variant_line.split('\t')
        allele_ref = fields[3]
        allele_alt = fields[4]
        parental_sup_info = fields[self.parental_sup_idx].split(':')
        parental_inf_info = fields[self.parental_inf_idx].split(':')
        pool_sup_info = fields[self.pool_sup_idx].split(':')
        pool_rnd_info = fields[self.pool_rnd_idx].split(':')
        sample_info = [parental_sup_info, parental_inf_info, pool_sup_info, pool_rnd_info]
        return allele_ref, allele_alt, sample_info

    def get_major_nuc(self, info: str) -> Tuple[str, List[int]]:
        '''Obtém o nucleotídeo com maior contagem'''
        nucs = ['A', 'C', 'G', 'T']
        counts = [int(x) for x in info.split(',')]
        i = counts.index(max(counts))
        return nucs[i], counts

    def get_counts(self, sample_info: List[str]) -> List[int]:
        '''Obtém contagens de um sample info'''
        counts = [int(x) for x in sample_info[4].split(',')]
        return counts

    def get_nuc_index(self, nuc: str) -> int:
        return ['A', 'C', 'G', 'T'].index(nuc)

    # Implementação dos filtros
    def filter_at_least(self, n_reads: int = 1) -> Tuple[List[str], List[str]]:
        '''Filtro: Contagem Mínima de Alelos de Referência'''
        output_lines = []
        filtered_lines = []
        n_lines = len(self.variant_lines)

        print(f'Filtering with at least {n_reads} reads')
        print(f'Number of variant lines: {n_lines}')

        for i, line in enumerate(self.variant_lines):
            allele_ref, allele_alt, sample_info = self.get_info_fields(line)
            print(f'{i} Allele Ref: {allele_ref}')
            print(f'{i} Allele Alt: {allele_alt}')
            print(f'{i} Sample Info: {sample_info}')

            parental_sup_info = sample_info[0]
            pool_sup_info = sample_info[2]

            print(f'{i} Parental Sup Info: {parental_sup_info}')
            print(f'{i} Pool Sup Info: {pool_sup_info}')

            print(f'Getting major nuc for parental sup info: {parental_sup_info[4]}')
            parental_sup_allele, parental_sup_counts = self.get_major_nuc(parental_sup_info[4])
            print(f'{i} Parental sup allele: {parental_sup_allele}')
            print(f'{i} Parental sup counts: {parental_sup_counts}')

            print(f'Getting major nuc for pool sup info: {pool_sup_info[4]}')
            pool_sup_allele, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])
            print(f'{i} Pool sup allele: {pool_sup_allele}')
            print(f'{i} Pool sup counts: {pool_sup_counts}')

            if parental_sup_allele != pool_sup_allele:
                filtered_lines.append(line)
                continue

            allele_freq_on_sup_pool = pool_sup_counts[self.get_nuc_index(parental_sup_allele)]
            if allele_freq_on_sup_pool >= n_reads:
                output_lines.append(line)
            else:
                filtered_lines.append(line)

        return output_lines, filtered_lines

    def filter_percent_threshold(self, threshold: float = 0.75) -> Tuple[List[str], List[str]]:
        '''Filtro: Limiar Percentual'''
        output_lines = []
        filtered_lines = []

        for line in self.variant_lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(line)
            parental_sup_info = sample_info[0]
            pool_sup_info = sample_info[2]

            pool_sup_counts = self.get_counts(pool_sup_info)
            pool_sup_total_reads = sum(pool_sup_counts)

            parental_sup_allele, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_allele, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])

            if parental_sup_allele != pool_sup_allele:
                filtered_lines.append(line)
                continue

            allele_freq_on_sup_pool = pool_sup_counts[self.get_nuc_index(parental_sup_allele)]
            pool_sup_allele_ratio = float(allele_freq_on_sup_pool) / float(pool_sup_total_reads)

            if pool_sup_allele_ratio >= threshold:
                output_lines.append(line)
            else:
                filtered_lines.append(line)

        return output_lines, filtered_lines

    def filter_parental_sup_greater(self) -> Tuple[List[str], List[str]]:
        '''Filtro: Alelo de Referência Dominante'''
        output_lines = []
        filtered_lines = []

        for line in self.variant_lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(line)
            parental_sup_info = sample_info[0]
            pool_sup_info = sample_info[2]

            parental_sup_allele, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_allele, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])

            if parental_sup_allele != pool_sup_allele:
                filtered_lines.append(line)
                continue

            if max(pool_sup_counts) == pool_sup_counts[self.get_nuc_index(pool_sup_allele)]:
                output_lines.append(line)
            else:
                filtered_lines.append(line)

        return output_lines, filtered_lines

    def filter_diff_from_greater(self, diff_max: float = 0.1) -> Tuple[List[str], List[str]]:
        '''Filtro: Diferença da Frequência Máxima'''
        output_lines = []
        filtered_lines = []

        for line in self.variant_lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(line)
            parental_sup_info = sample_info[0]
            pool_sup_info = sample_info[2]

            parental_sup_allele, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_allele, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])
            pool_sup_total_reads = sum(pool_sup_counts)

            pool_sup_percents = [float(x) / float(pool_sup_total_reads) for x in pool_sup_counts]
            greater_perc = max(pool_sup_percents)
            pool_sup_allele_perc = pool_sup_percents[self.get_nuc_index(parental_sup_allele)]
            diff = abs(greater_perc - pool_sup_allele_perc)

            if diff_max >= diff >= 0.0:
                output_lines.append(line)
            else:
                filtered_lines.append(line)

        return output_lines, filtered_lines

    def filter_rnd_mean(self, threshold: float = 0.75, avg: float = 0.5, std: float = 0.02) -> Tuple[
        List[str], List[str]]:
        '''Filtro: Média do Pool Aleatório'''
        output_lines = []
        filtered_lines = []

        for line in self.variant_lines:
            allele_ref, allele_alt, sample_info = self.get_info_fields(line)
            parental_sup_info = sample_info[0]
            pool_sup_info = sample_info[2]
            pool_rnd_info = sample_info[3]

            parental_sup_allele, _ = self.get_major_nuc(parental_sup_info[4])
            pool_sup_allele, pool_sup_counts = self.get_major_nuc(pool_sup_info[4])
            pool_rnd_counts = self.get_counts(pool_rnd_info)

            if parental_sup_allele != pool_sup_allele:
                filtered_lines.append(line)
                continue

            pool_sup_total_reads = sum(pool_sup_counts)
            pool_rnd_total_reads = sum(pool_rnd_counts)
            pool_sup_percents = [float(x) / float(pool_sup_total_reads) for x in pool_sup_counts]
            pool_rnd_percents = [float(x) / float(pool_rnd_total_reads) for x in pool_rnd_counts]

            ref_sup_perc = pool_sup_percents[self.get_nuc_index(parental_sup_allele)]
            ref_rnd_perc = pool_rnd_percents[self.get_nuc_index(parental_sup_allele)]

            cond1 = ref_sup_perc >= threshold
            cond2 = (avg + std) >= ref_rnd_perc >= (avg - std)

            if cond1 and cond2:
                output_lines.append(line)
            else:
                filtered_lines.append(line)

        return output_lines, filtered_lines


class VCFProcessorGUI:
    def __init__(self, root, logger=None):
        self.root = root
        self.logger = logger
        self.root.title("Processador de VCF")
        self.root.geometry("900x800")

        self.processor = None
        self.setup_logging()
        self.setup_ui()

    def setup_logging(self):
        '''Configura o sistema de logging'''
        self.logger = logging.getLogger('VCFProcessor')
        self.logger.setLevel(logging.INFO)

    def setup_ui(self):
        # Frame principal
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configurar grid
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # Seção de arquivo VCF
        ttk.Label(main_frame, text="Arquivo VCF:", font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky=tk.W,
                                                                                    pady=5)

        file_frame = ttk.Frame(main_frame)
        file_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        file_frame.columnconfigure(0, weight=1)

        self.file_entry = ttk.Entry(file_frame)
        self.file_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(file_frame, text="Procurar", command=self.browse_file).grid(row=0, column=1)

        # Botão para carregar VCF
        ttk.Button(main_frame, text="Carregar VCF", command=self.load_vcf).grid(row=2, column=0, columnspan=2, pady=5)

        # Seção de configuração das colunas
        self.column_config_frame = ttk.LabelFrame(main_frame, text="Configuração das Colunas de Amostras", padding="10")
        self.column_config_frame.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        self.column_config_frame.columnconfigure(1, weight=1)

        # Inicialmente escondido até carregar o VCF
        self.column_config_frame.grid_remove()

        # Variáveis para os Comboboxes
        self.parental_sup_var = tk.StringVar()
        self.parental_inf_var = tk.StringVar()
        self.pool_sup_var = tk.StringVar()
        self.pool_rnd_var = tk.StringVar()

        # Labels e Comboboxes para configuração das colunas
        ttk.Label(self.column_config_frame, text="Parental Superior:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.parental_sup_combo = ttk.Combobox(self.column_config_frame, textvariable=self.parental_sup_var,
                                               state="readonly")
        self.parental_sup_combo.grid(row=0, column=1, sticky=(tk.W, tk.E), pady=5, padx=(10, 0))

        ttk.Label(self.column_config_frame, text="Parental Inferior:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.parental_inf_combo = ttk.Combobox(self.column_config_frame, textvariable=self.parental_inf_var,
                                               state="readonly")
        self.parental_inf_combo.grid(row=1, column=1, sticky=(tk.W, tk.E), pady=5, padx=(10, 0))

        ttk.Label(self.column_config_frame, text="Pool Superior:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.pool_sup_combo = ttk.Combobox(self.column_config_frame, textvariable=self.pool_sup_var, state="readonly")
        self.pool_sup_combo.grid(row=2, column=1, sticky=(tk.W, tk.E), pady=5, padx=(10, 0))

        ttk.Label(self.column_config_frame, text="Pool Aleatório:").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.pool_rnd_combo = ttk.Combobox(self.column_config_frame, textvariable=self.pool_rnd_var, state="readonly")
        self.pool_rnd_combo.grid(row=3, column=1, sticky=(tk.W, tk.E), pady=5, padx=(10, 0))

        # Botão para confirmar configuração
        ttk.Button(self.column_config_frame, text="Confirmar Configuração", command=self.confirm_column_config).grid(
            row=4, column=0, columnspan=2, pady=10)

        # Seção de seleção de filtro
        ttk.Label(main_frame, text="Filtro:", font=('Arial', 10, 'bold')).grid(row=4, column=0, sticky=tk.W, pady=10)

        self.filter_var = tk.StringVar()
        self.filter_combo = ttk.Combobox(main_frame, textvariable=self.filter_var, state="readonly")
        self.filter_combo['values'] = (
            "Contagem Mínima de Alelos",
            "Limiar Percentual",
            "Alelo de Referência Dominante",
            "Diferença da Frequência Máxima",
            "Média do Pool Aleatório"
        )
        self.filter_combo.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        self.filter_combo.bind('<<ComboboxSelected>>', self.on_filter_select)

        # Frame de parâmetros
        self.param_frame = ttk.LabelFrame(main_frame, text="Parâmetros do Filtro", padding="10")
        self.param_frame.grid(row=6, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        self.param_frame.columnconfigure(1, weight=1)

        # Área de logs
        log_frame = ttk.LabelFrame(main_frame, text="Logs", padding="5")
        log_frame.grid(row=7, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        main_frame.rowconfigure(7, weight=1)

        self.log_text = tk.Text(log_frame, height=15, wrap=tk.WORD)
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)

        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # Frame de botões
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=8, column=0, columnspan=2, pady=10)

        ttk.Button(button_frame, text="Aplicar Filtro", command=self.apply_filter).pack(side=tk.LEFT, padx=5)
        self.save_button = ttk.Button(button_frame, text="Salvar Resultado", command=self.save_result, state="disabled")
        self.save_button.pack(side=tk.LEFT, padx=5)

        self.current_output_lines = []
        self.log("Aplicação iniciada. Selecione um arquivo VCF.")

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Selecionar arquivo VCF",
            filetypes=[("VCF files", "*.vcf"), ("All files", "*.*")]
        )
        if filename:
            self.file_entry.delete(0, tk.END)
            self.file_entry.insert(0, filename)
            self.log(f"Arquivo selecionado: {os.path.basename(filename)}")

    def load_vcf(self):
        file_path = self.file_entry.get()
        if not file_path or not os.path.exists(file_path):
            messagebox.showerror("Erro", "Por favor, selecione um arquivo VCF válido.")
            return

        try:
            print(f'Initializing VCFProcessor with logger: {self.logger}')
            self.processor = VCFProcessor(self)
            self.processor.set_file_path(file_path)

            if not self.processor.load_vcf_file():
                return

            # Mostrar frame de configuração de colunas
            self.column_config_frame.grid()

            # Preencher comboboxes com as colunas detectadas
            sample_columns = self.processor.sample_columns
            self.parental_sup_combo['values'] = sample_columns
            self.parental_inf_combo['values'] = sample_columns
            self.pool_sup_combo['values'] = sample_columns
            self.pool_rnd_combo['values'] = sample_columns

            # Selecionar valores padrão (primeiras 4 colunas na ordem)
            if len(sample_columns) >= 4:
                self.parental_sup_var.set(sample_columns[0])
                self.parental_inf_var.set(sample_columns[1])
                self.pool_sup_var.set(sample_columns[2])
                self.pool_rnd_var.set(sample_columns[3])

            # self.log("Arquivo VCF carregado com sucesso!")
            self.log(f"Colunas de amostras detectadas: {', '.join(sample_columns)}")
            self.log("Por favor, configure as colunas acima e clique em 'Confirmar Configuração'")

        except Exception as e:
            self.log(f"Erro ao carregar VCF: {str(e)}")
            messagebox.showerror("Erro", f"Erro ao carregar VCF: {str(e)}")

    def confirm_column_config(self):
        '''Confirma a configuração das colunas selecionadas pelo usuário'''
        if not self.processor:
            messagebox.showerror("Erro", "Nenhum arquivo VCF carregado.")
            return

        parental_sup = self.parental_sup_var.get()
        parental_inf = self.parental_inf_var.get()
        pool_sup = self.pool_sup_var.get()
        pool_rnd = self.pool_rnd_var.get()

        # Verificar se todas as colunas foram selecionadas e são únicas
        selected_columns = [parental_sup, parental_inf, pool_sup, pool_rnd]
        if not all(selected_columns):
            messagebox.showerror("Erro", "Por favor, selecione todas as colunas.")
            return

        if len(set(selected_columns)) != 4:
            messagebox.showerror("Erro", "Por favor, selecione colunas diferentes para cada tipo.")
            return

        try:
            if self.processor.set_sample_columns(parental_sup, parental_inf, pool_sup, pool_rnd):
                self.log("Configuração de colunas confirmada com sucesso!")
                self.log(f"Parental Superior: {parental_sup}")
                self.log(f"Parental Inferior: {parental_inf}")
                self.log(f"Pool Superior: {pool_sup}")
                self.log(f"Pool Aleatório: {pool_rnd}")

                # Realizar a correção dos genótipos
                self.processor.correct_genotypes(verbose=True)
            else:
                messagebox.showerror("Erro", "Erro ao configurar colunas.")

        except Exception as e:
            self.log(f"Erro ao configurar colunas: {str(e)}")
            messagebox.showerror("Erro", f"Erro ao configurar colunas: {str(e)}")

    def on_filter_select(self, event):
        selected_filter = self.filter_combo.get()
        self.clear_param_frame()

        if selected_filter == "Contagem Mínima de Alelos":
            self.create_at_least_params()
        elif selected_filter == "Limiar Percentual":
            self.create_percent_threshold_params()
        elif selected_filter == "Diferença da Frequência Máxima":
            self.create_diff_from_greater_params()
        elif selected_filter == "Média do Pool Aleatório":
            self.create_rnd_mean_params()

    def clear_param_frame(self):
        for widget in self.param_frame.winfo_children():
            widget.destroy()

    def create_at_least_params(self):
        ttk.Label(self.param_frame, text="Contagem mínima de leituras:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.n_reads_var = tk.StringVar(value="1")
        ttk.Entry(self.param_frame, textvariable=self.n_reads_var, width=10).grid(row=0, column=1, sticky=tk.W, pady=2)

    def create_percent_threshold_params(self):
        ttk.Label(self.param_frame, text="Limiar percentual (%):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.threshold_var = tk.StringVar(value="75")
        ttk.Entry(self.param_frame, textvariable=self.threshold_var, width=10).grid(row=0, column=1, sticky=tk.W,
                                                                                    pady=2)

    def create_diff_from_greater_params(self):
        ttk.Label(self.param_frame, text="Diferença máxima (%):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.diff_max_var = tk.StringVar(value="10")
        ttk.Entry(self.param_frame, textvariable=self.diff_max_var, width=10).grid(row=0, column=1, sticky=tk.W, pady=2)

    def create_rnd_mean_params(self):
        ttk.Label(self.param_frame, text="Limiar (%):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.rnd_threshold_var = tk.StringVar(value="75")
        ttk.Entry(self.param_frame, textvariable=self.rnd_threshold_var, width=10).grid(row=0, column=1, sticky=tk.W,
                                                                                        pady=2)

        ttk.Label(self.param_frame, text="Média (%):").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.avg_count_var = tk.StringVar(value="50")
        ttk.Entry(self.param_frame, textvariable=self.avg_count_var, width=10).grid(row=1, column=1, sticky=tk.W,
                                                                                    pady=2)

        ttk.Label(self.param_frame, text="Desvio padrão (%):").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.std_dev_var = tk.StringVar(value="2")
        ttk.Entry(self.param_frame, textvariable=self.std_dev_var, width=10).grid(row=2, column=1, sticky=tk.W, pady=2)

    def log(self, message: str):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def apply_filter(self):
        if not self.processor or not self.processor.variant_lines:
            messagebox.showerror("Erro", "Por favor, carregue e configure um arquivo VCF primeiro.")
            return

        if not self.filter_var.get():
            messagebox.showerror("Erro", "Por favor, selecione um filtro.")
            return

        # Verificar se as colunas foram configuradas
        if any(idx is None for idx in [self.processor.parental_sup_idx, self.processor.parental_inf_idx,
                                       self.processor.pool_sup_idx, self.processor.pool_rnd_idx]):
            messagebox.showerror("Erro", "Por favor, configure as colunas primeiro.")
            return

        try:
            print(f'Filter selected: {self.filter_var.get()}')
            selected_filter = self.filter_combo.get()
            original_count = len(self.processor.variant_lines)

            self.log(f"Aplicando filtro: {selected_filter}")
            self.log(f"Número original de variantes: {original_count}")

            if selected_filter == "Contagem Mínima de Alelos":
                n_reads = int(self.n_reads_var.get())
                output_lines, filtered_lines = self.processor.filter_at_least(n_reads)

            elif selected_filter == "Limiar Percentual":
                threshold = float(self.threshold_var.get()) / 100.0
                output_lines, filtered_lines = self.processor.filter_percent_threshold(threshold)

            elif selected_filter == "Alelo de Referência Dominante":
                output_lines, filtered_lines = self.processor.filter_parental_sup_greater()

            elif selected_filter == "Diferença da Frequência Máxima":
                diff_max = float(self.diff_max_var.get()) / 100.0
                output_lines, filtered_lines = self.processor.filter_diff_from_greater(diff_max)

            elif selected_filter == "Média do Pool Aleatório":
                threshold = float(self.rnd_threshold_var.get()) / 100.0
                avg = float(self.avg_count_var.get()) / 100.0
                std = float(self.std_dev_var.get()) / 100.0
                output_lines, filtered_lines = self.processor.filter_rnd_mean(threshold, avg, std)

            filtered_count = len(output_lines)
            removed_count = original_count - filtered_count

            self.log("Filtro aplicado com sucesso!")
            self.log(f"Variantes originais: {original_count}")
            self.log(f"Variantes após filtro: {filtered_count}")
            self.log(f"Variantes removidas: {removed_count}")
            self.log(f"Taxa de retenção: {filtered_count / original_count * 100:.2f}%")

            self.current_output_lines = output_lines
            self.save_button.config(state="normal")

        except Exception as e:
            self.log(f"Erro ao aplicar filtro: {str(e)}")
            messagebox.showerror("Erro", f"Erro ao aplicar filtro: {str(e)}")

    def save_result(self):
        if not self.current_output_lines or not self.processor:
            messagebox.showwarning("Aviso", "Nenhum resultado para salvar.")
            return

        default_name = f"filtered_{os.path.basename(self.processor.vcf_file_path)}"
        filename = filedialog.asksaveasfilename(
            title="Salvar arquivo VCF filtrado",
            initialfile=default_name,
            defaultextension=".vcf",
            filetypes=[("VCF files", "*.vcf"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w') as f:
                    # Escrever headers
                    for header in self.processor.metadata_lines:
                        f.write(header)
                    # Escrever variantes filtradas
                    for line in self.current_output_lines:
                        f.write(line)

                self.log(f"Arquivo salvo: {filename}")
                messagebox.showinfo("Sucesso",
                                    f"Arquivo salvo com sucesso!\n{len(self.current_output_lines)} variantes exportadas.")

            except Exception as e:
                self.log(f"Erro ao salvar arquivo: {str(e)}")
                messagebox.showerror("Erro", f"Erro ao salvar arquivo: {str(e)}")


def main():
    root = tk.Tk()
    app = VCFProcessorGUI(root, logging)
    root.mainloop()


if __name__ == "__main__":
    main()
