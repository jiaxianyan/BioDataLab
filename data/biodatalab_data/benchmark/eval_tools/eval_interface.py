from .common_eval import  *  
from .specific_eval import *

eval_tool_dict = {
                    'eval_str_from_file_equal':eval_str_from_file_equal,
                    'eval_fasta_file_equal': eval_fasta_from_file_equal,

                    'eval_numeric_value_from_file_equal_round': eval_numeric_value_from_file_equal_round,
                    'eval_numeric_value_list_from_file_equal_round': eval_numeric_value_list_from_file_equal_round,

                    'eval_str_list_from_json_equal_unsort': eval_str_list_from_json_equal_unsort,
                    'eval_str_list_from_json_equal_sort': eval_str_list_from_json_equal_sort,

                    'eval_dict_from_json_equal': eval_dict_from_json_equal,

                    'eval_csv_file_equal_sort': eval_csv_file_equal_sort,
                    'eval_csv_file_equal_unsort': eval_csv_file_equal_unsort,

                    'eval_sam_bam_equal': eval_sam_bam_flagstat,

                    'eval_tsv_file_unsort': eval_tsv_file_unsort,
                    'eval_tsv_file_sort': eval_tsv_file_sort,

                    'eval_bed_file_equal': eval_bed_file_equal,

                    'eval_vcf_file_equal': eval_vcf_file_equal,

                    'eval_fastq_file_equal': eval_fastq_file_equal,

                    'eval_vcf_gz_file_equal': eval_vcf_gz_file_equal,

                    'eval_gtf_file_equal': eval_gtf_file_equal,

                    'eval_fq_gz_file_equal': eval_fq_gz_file_equal,

                    'eval_psl_file_equal': eval_psl_file_equal,

                    'eval_rds_file_equal': eval_rds_file_equal,

                    'eval_h5ad_file_equal': eval_h5ad_file_equal,

                    'eval_specific_coding_10': eval_specific_coding_10
                    }