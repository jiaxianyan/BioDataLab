# cmw
from .adcdb_extract_2 import adcdb_extract_2
from .asmdb_retrieval import asmdb_retrieval
from .bioka_retrieval import bioka_retrieval
from .cds_db_retrieval import cds_db_retrieval
from .circmine_retrieval import circmine_retrieval
from .crost_retrieval import crost_retrieval
from .ctr_db_retrieval import ctr_db_retrieval
from .ddinter_retrieval import ddinter_retrieval
from .fusionneoantigen_extract import fusionneoantigen_extract
from .npcdr_extract_1 import npcdr_extract_1
from .npcdr_retrieval import npcdr_retrieval
from .ravar_extract import ravar_extract
from .scan_retrieval import scan_retrieval
from .scqtlbase_retrieval import scqtlbase_retrieval
from .stemdriver_retrieval import stemdriver_retrieval
from .tf_marker_retrieval import tf_marker_retrieval
from .themarker_extract_1 import themarker_extract_1
from .adcdb_extract_1 import adcdb_extract_1
from .circmine_extract import circmine_extract
from .cyanoomicsdb_retrieval_2 import cyanoomicsdb_retrieval_2
from .dda_extract import dda_extract
from .disco_extract import disco_extract
from .dntppooldb_extract import dntppooldb_extract
from .macc_extract import macc_extract
from .npcdr_extract_2 import npcdr_extract_2
from .pcmdb_extract import pcmdb_extract
from .pharmgwas_extract import pharmgwas_extract
from .pronab_extract import pronab_extract
from .themarker_extract_2 import themarker_extract_2

# jxy
from .covpdb_retrieval import covpdb_retrieval
from .covpdb_integration import covpdb_integration
from .covpdb_annotate import covpdb_annotate
from .scov2_md_annotate import scov2_md_annotate
from .scovid_refinement import scovid_refinement
from .vareps_annotate import vareps_annotation
from .cyanoomicsdb_annotate_1 import cyanoomicsdb_annotate_1
from .mvip_annotate import mvip_annotate
from .pncshub_annotate import pncshub_annotate
from .compodynamics_integration import compodynamics_integration
from .compodynamics_annotate import compodynamics_annotate
from .metazexp_refinement import metazexp_refinement  
from .metazexp_annotate import metazexp_annotate
from .diana_mited_refinement import diana_mited_refinement
from .cellstar_integration import cellstar_integration
from .colocdb_refinement import colocdb_refinement
from .pgs_depot_refinement import pgs_depot_refinement
from .ravar_refinement_1 import ravar_refinement_1
from .ravar_refinement_2 import ravar_refinement_2
from .scqtlbase_refinement import scqtlbase_refinement
from .themarker_annotate import themarker_annotate
from .cancerscem_annotate import cancerscem_annotate
from .ddinter_integration_1 import ddinter_integration_1
from .ddinter_integration_2 import ddinter_integration_2
from .ddinter_annotate_1 import ddinter_annotate_1
from .ddinter_annotate_2 import ddinter_annotate_2
from .tf_marker_annotate import tf_marker_annotate
from .plantpad_annotate import plantpad_annotate
from .covid_19_extract import covid_19_extract
from .amdb_retrieval import amdb_retrieval
from .amdb_extract import amdb_extract
from .vimic_extract import vimic_extract
from .zover_extract import zover_extract
from .pcmdb_extract_2 import pcmdb_extract_2
from .covid_19_integration import covid_19_integration
from .cancerscem_annotate_2 import cancerscem_annotate_2
from .dda_refinement import dda_refinement
from .fusionneoantigen_annotate import fusionneoantigen_annotate
from .disco_refinement import disco_refinement
from .gpedit_refinement import gpedit_refinement
from .dntppooldb_refinement import dntppooldb_refinement
from .clinicalomicsdb_annotate import clinicalomicsdb_annotate
from .oncodb_annotate import oncodb_annotate
from .cancerproteome_annotate import cancerproteome_annotate
from .cellcommunet_refinement import cellcommunet_refinement
from .cancerproteome_annotate_2 import cancerproteome_annotate_2
from .scapaatlas_annotate import scapaatlas_annotate
from .a3d_modb_retrieval import a3d_modb_retrieval
from .cyanoomicsdb_annotate_2 import cyanoomicsdb_annotate_2
from .m2or_annotate import m2or_annotate
from .cancermirnome_annotate import cancermirnome_annotate
from .drmref_annotate import drmref_annotate
from .diana_mited_retrieval import diana_mited_retrieval
from .bioka_extract import bioka_extract
from .asmdb_refinement_1 import asmdb_refinement_1
from .asmdb_refinement_2 import asmdb_refinement_2
from .asmdb_annotate import asmdb_annotate
from .m2or_refinement import m2or_refinement
from .mbodymap_integration import mbodymap_integration
from .fusionneoantigen_annotate_2 import fusionneoantigen_annotate_2

# zjt 
from .plabdab_retrieval import plabdab_retrieval
from .plabdab_annotate_1 import plabdab_annotate_1
from .plabdab_annotate_2 import plabdab_annotate_2
from .atlas_retrieval_1 import atlas_retrieval_1
from .atlas_retrieval_2 import atlas_retrieval_2
from .inclusive_retrieval import inclusive_retrieval
from .inclusive_extract_1 import inclusive_extract_1
from .inclusive_extract_2 import inclusive_extract_2
from .inclusive_extract_3 import inclusive_extract_3
from .kincore_retrieval import kincore_retrieval
from .kincore_renumbering import kincore_renumbering


eval_tool_dict = {
                    'adcdb_extract_2': adcdb_extract_2,
                    'asmdb_retrieval': asmdb_retrieval,
                    'bioka_retrieval': bioka_retrieval,
                    'cds_db_retrieval': cds_db_retrieval,
                    'circmine_retrieval': circmine_retrieval,
                    'crost_retrieval': crost_retrieval,
                    'ctr_db_retrieval': ctr_db_retrieval,
                    'ddinter_retrieval': ddinter_retrieval,
                    'fusionneoantigen_extract': fusionneoantigen_extract,
                    'npcdr_extract_1': npcdr_extract_1,
                    'npcdr_retrieval': npcdr_retrieval,
                    'ravar_extract': ravar_extract,
                    'scan_retrieval': scan_retrieval,
                    'scqtlbase_retrieval': scqtlbase_retrieval,
                    'stemdriver_retrieval': stemdriver_retrieval,
                    'tf_marker_retrieval': tf_marker_retrieval,
                    'themarker_extract_1': themarker_extract_1,
                    'adcdb_extract_1': adcdb_extract_1,
                    'circmine_extract': circmine_extract,
                    'cyanoomicsdb_retrieval_2': cyanoomicsdb_retrieval_2,
                    'dda_extract': dda_extract,
                    'disco_extract': disco_extract,
                    'dntppooldb_extract': dntppooldb_extract,
                    'macc_extract': macc_extract,
                    'npcdr_extract_2': npcdr_extract_2,
                    'pcmdb_extract': pcmdb_extract,
                    'pharmgwas_extract': pharmgwas_extract,
                    'pronab_extract': pronab_extract,
                    'themarker_extract_2': themarker_extract_2,
                    'dntppooldb_refinement': dntppooldb_refinement,
                    
                    'covpdb_retrieval': covpdb_retrieval,
                    'covpdb_integration': covpdb_integration,
                    'covpdb_annotate': covpdb_annotate,
                    'scov2_md_annotate': scov2_md_annotate,
                    'scovid_refinement': scovid_refinement,
                    'vareps_annotation': vareps_annotation,
                    'cyanoomicsdb_annotate_1': cyanoomicsdb_annotate_1,
                    'mvip_annotate': mvip_annotate,
                    'pncshub_annotate': pncshub_annotate,
                    'compodynamics_integration': compodynamics_integration,
                    'compodynamics_annotate': compodynamics_annotate,
                    'metazexp_refinement': metazexp_refinement,
                    'metazexp_annotate': metazexp_annotate,
                    'diana_mited_refinement': diana_mited_refinement,
                    'cellstar_integration': cellstar_integration,
                    'colocdb_refinement': colocdb_refinement,
                    'pgs_depot_refinement': pgs_depot_refinement,
                    'ravar_refinement_1': ravar_refinement_1,
                    'ravar_refinement_2': ravar_refinement_2,
                    'scqtlbase_refinement': scqtlbase_refinement,
                    'themarker_annotate': themarker_annotate,
                    'cancerscem_annotate': cancerscem_annotate,
                    'ddinter_integration_1': ddinter_integration_1,
                    'ddinter_integration_2': ddinter_integration_2,
                    'ddinter_annotate_1': ddinter_annotate_1,
                    'ddinter_annotate_2': ddinter_annotate_2,
                    'tf_marker_annotate': tf_marker_annotate,
                    'plantpad_annotate': plantpad_annotate,
                    'covid_19_extract': covid_19_extract,
                    'amdb_retrieval': amdb_retrieval,
                    'amdb_extract': amdb_extract,
                    'vimic_extract': vimic_extract,
                    'zover_extract': zover_extract,
                    'pcmdb_extract_2': pcmdb_extract_2,
                    'covid_19_integration': covid_19_integration,
                    'cancerscem_annotate_2': cancerscem_annotate_2,
                    'dda_refinement': dda_refinement,
                    'fusionneoantigen_annotate': fusionneoantigen_annotate,
                    'disco_refinement': disco_refinement,
                    'gpedit_refinement': gpedit_refinement,
                    'clinicalomicsdb_annotate': clinicalomicsdb_annotate,
                    'oncodb_annotate': oncodb_annotate,
                    'cancerproteome_annotate': cancerproteome_annotate,
                    'cellcommunet_refinement': cellcommunet_refinement,
                    'cancerproteome_annotate_2': cancerproteome_annotate_2,
                    'scapaatlas_annotate': scapaatlas_annotate,
                    'a3d_modb_retrieval': a3d_modb_retrieval,
                    'cyanoomicsdb_annotate_2': cyanoomicsdb_annotate_2,
                    'm2or_annotate': m2or_annotate,
                    'cancermirnome_annotate': cancermirnome_annotate,
                    'drmref_annotate': drmref_annotate,
                    'diana_mited_retrieval': diana_mited_retrieval,
                    'bioka_extract': bioka_extract,
                    'asmdb_refinement_1': asmdb_refinement_1,
                    'asmdb_refinement_2': asmdb_refinement_2,
                    'asmdb_annotate': asmdb_annotate,
                    'm2or_refinement': m2or_refinement,
                    'mbodymap_integration': mbodymap_integration,
                    'fusionneoantigen_annotate_2': fusionneoantigen_annotate_2,
                    
                    # zjt works
                    'plabdab_retrieval': plabdab_retrieval,
                    'plabdab_annotate_1': plabdab_annotate_1,
                    'plabdab_annotate_2': plabdab_annotate_2,
                    'atlas_retrieval_1': atlas_retrieval_1,
                    'atlas_retrieval_2': atlas_retrieval_2,
                    'inclusive_retrieval': inclusive_retrieval,
                    'inclusive_extract_1': inclusive_extract_1,
                    'inclusive_extract_2': inclusive_extract_2,
                    'inclusive_extract_3': inclusive_extract_3,
                    'kincore_retrieval': kincore_retrieval,
                    'kincore_renumbering': kincore_renumbering,
                                        
                 }

