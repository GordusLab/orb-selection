library(here)

source(here("scripts/05_permulation_loss_dup/permulations.R"))

# # Run excluding araneids

# testCatPerms <- run_categorical_permulations(
#         foreground_list_filename = here("data/orbweavers-list.txt"),
#         ntrees = 10000,
#         save_rdata_path = here("data/perms10000_no_araneids.RData"),
#         excluded_tips = c(
#             "Drosophila_melanogaster",
#             "Antrodiaetus_roretzi",
#             "Orchestina_okitsui",
#             "Falcileptoneta_japonica",
#             "Masirana_silvicola",
#             "Phonognatha_graeffei",
#             "Caerostris_extrusa",
#             "Caerostris_darwini",
#             "Araneus_inustus",
#             "Cyrtophora_cicatrosa",
#             "Argiope_bruennichi",
#             "Argiope_keyserlingi",
#             "Argiope_aetheroides",
#             "Metazygia_zilloides",
#             "Acrosomoides_sp_IDV7426",
#             "Backobourkia_brouni",
#             "Dolophones_sp_IDV6683",
#             "Cyclosa_octotuberculata",
#             "Plebs_sachalinensis",
#             "Acroaspis_sp_IDV6688",
#             "Eriophora_pustulosa",
#             "Neoscona_subpullata",
#             "Neoscona_scylla",
#             "Neoscona_theisi",
#             "Neoscona_adianta",
#             "Neoscona_nautica",
#             "Poecilopachys_australasia",
#             "Hypsosinga_pygmaea",
#             "Lariniaria_argiopiformis",
#             "Larinia_phthisica",
#             "Nuctenea_umbratica",
#             "Larinioides_cornutus", 
#             "Ordgarius_hobsoni",
#             "Ordgarius_sexspinosus",
#             "Cyrtarachne_bufo"
#         )
# )

# save_tip_values(testCatPerms, here("data/perms_tip_values_no_araneids.csv"))

# # Run excluding tetragnathids

# testCatPerms <- run_categorical_permulations(
#         foreground_list_filename = here("data/orbweavers-list.txt"),
#         ntrees = 10000,
#         save_rdata_path = here("data/perms10000_no_tetragnathids.RData"),
#         excluded_tips = c(
#             "Drosophila_melanogaster",
#             "Antrodiaetus_roretzi",
#             "Orchestina_okitsui",
#             "Falcileptoneta_japonica",
#             "Masirana_silvicola",
#             "Metleucauge_yunohamensis",
#             "Mesida_sp_IDV5268",
#             "Leucauge_dromedari",
#             "Tetragnatha_yesoensis",
#             "Tetragnatha_vermiformis",
#             "Tetragnatha_maxillosa",
#             "Tetragnatha_mandibulata",
#             "Tetragnatha_extensa",
#             "Tetragnatha_quasimodo",
#             "Tetragnatha_brevignatha",
#             "Tetragnatha_waikamoi",
#             "Tetragnatha_filiciphilia",
#             "Tetragnatha_stelarobusta",
#             "Tetragnatha_paludicola"
#         )
# )

# save_tip_values(testCatPerms, here("data/perms_tip_values_no_tetragnathids.csv"))

# Run only araneids
# For odds ratio test, use "non-orb-weavers-list.txt" as background and "araneids-orbweavers.txt" as foreground
testCatPerms <- run_categorical_permulations(
        foreground_list_filename = here("data/orbweavers-list.txt"),
        ntrees = 10000,
        save_rdata_path = here("data/perms10000_only_araneids.RData"),
        excluded_tips = c(
            "Drosophila_melanogaster",
            "Antrodiaetus_roretzi",
            "Orchestina_okitsui",
            "Falcileptoneta_japonica",
            "Masirana_silvicola",
            "Metleucauge_yunohamensis",
            "Mesida_sp_IDV5268",
            "Leucauge_dromedari",
            "Tetragnatha_yesoensis",
            "Tetragnatha_vermiformis",
            "Tetragnatha_maxillosa",
            "Tetragnatha_mandibulata",
            "Tetragnatha_extensa",
            "Tetragnatha_quasimodo",
            "Tetragnatha_brevignatha",
            "Tetragnatha_waikamoi",
            "Tetragnatha_filiciphilia",
            "Tetragnatha_stelarobusta",
            "Tetragnatha_paludicola",
            "Uloborus_diversus",
            "Deinopis_sp_IDV7365",
            "Deinopis_sp_IDV6783"
        )
)

save_tip_values(testCatPerms, here("data/perms_tip_values_only_araneids.csv"))

# Run only tetragnathids
# For odds ratio test, use "non-orb-weavers-no-araneids.txt" as background and "tetragnathids.txt" as foreground
testCatPerms <- run_categorical_permulations(
        foreground_list_filename = here("data/orbweavers-list.txt"),
        ntrees = 10000,
        save_rdata_path = here("data/perms10000_only_tetragnathids.RData"),
        excluded_tips = c(
            "Drosophila_melanogaster",
            "Antrodiaetus_roretzi",
            "Orchestina_okitsui",
            "Falcileptoneta_japonica",
            "Masirana_silvicola",
            "Uloborus_diversus",
            "Deinopis_sp_IDV7365",
            "Deinopis_sp_IDV6783",
            "Phonognatha_graeffei",
            "Caerostris_extrusa",
            "Caerostris_darwini",
            "Araneus_inustus",
            "Cyrtophora_cicatrosa",
            "Argiope_bruennichi",
            "Argiope_keyserlingi",
            "Argiope_aetheroides",
            "Metazygia_zilloides",
            "Acrosomoides_sp_IDV7426",
            "Backobourkia_brouni",
            "Dolophones_sp_IDV6683",
            "Cyclosa_octotuberculata",
            "Plebs_sachalinensis",
            "Acroaspis_sp_IDV6688",
            "Eriophora_pustulosa",
            "Neoscona_subpullata",
            "Neoscona_scylla",
            "Neoscona_theisi",
            "Neoscona_adianta",
            "Neoscona_nautica",
            "Poecilopachys_australasia",
            "Hypsosinga_pygmaea",
            "Lariniaria_argiopiformis",
            "Larinia_phthisica",
            "Nuctenea_umbratica",
            "Larinioides_cornutus", 
            "Ordgarius_hobsoni",
            "Ordgarius_sexspinosus",
            "Cyrtarachne_bufo"
        )
)

save_tip_values(testCatPerms, here("data/perms_tip_values_only_tetragnathids.csv"))

