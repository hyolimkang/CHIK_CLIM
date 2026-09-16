# Pernambuco spatial model -- Section 1: official spatial partition.
#
# Builds PE_spatial_stratum_lookup.csv assigning every PE municipality to
# (a) its official health macroregion and (b) the final 5-way model stratum
# used in Spatial Model 1. The partition is administrative, defined
# independently of observed chikungunya incidence/wave timing/hotspot
# status -- per instruction, NOT derived from the spatial-turnover
# diagnostic.
#
# SOURCE (authoritative, official, external to this repo -- no PE
# macroregion lookup previously existed in-repo; only Ceara's did):
#   1. GERES (Gerencia Regional de Saude) -> municipality lists: transcribed
#      from "Enderecos e Municipios Geres PE-21.05.24.pdf", published by
#      APEVISA (Agencia Pernambucana de Vigilancia Sanitaria, an SES-PE
#      linked state agency), https://apevisa.pe.gov.br/ -- cross-checked
#      against the same list on SES-PE's own
#      https://cartadeservicos.saude.pe.gov.br/endereco-geres/ (accessed via
#      WebFetch 2026-09-14; both government sources agree verbatim).
#   2. GERES -> official macroregion grouping (Metropolitana = I+II+III+XII;
#      Agreste = IV+V; Sertao = VI+X+XI; Vale do Sao Francisco e Araripe =
#      VII+VIII+IX), taken verbatim from Pernambuco's official Plano
#      Estadual de Saude (PES) 2016-2019, published by CONASS
#      (https://www.conass.org.br/pdf/planos-estaduais-de-saude/PE_PES-2016-2019-FINAL_23_12_2016-1.pdf),
#      p.36 and section 2.5.3, which quotes the state's 2011 Plano Diretor
#      de Regionalizacao (PDR) revision verbatim ("Em 2011, foi realizada a
#      revisao do Plano Diretor de Regionalizacao (PDR) que definiu as
#      seguintes macrorregioes: Metropolitana: I, II, III e XII Regioes de
#      Saude [...]; Agreste: IV e V Regioes de Saude [...]; Sertao: VI, X e
#      XI Regioes de Saude [...]; Vale do Sao Francisco: VII, VIII e IX
#      Regioes De Saude [...]").
# This partition is NOT derived from observed chikungunya data.

required_packages <- c("dplyr", "readr", "stringr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Official GERES -> municipality lists (verbatim, ASCII-normalised for matching) ----
geres_list <- list(
  "I"   = c("Abreu e Lima","Aracoiaba","Cabo de Santo Agostinho","Camaragibe","Cha Grande",
            "Cha de Alegria","Gloria de Goita","Fernando de Noronha","Igarassu","Ipojuca",
            "Itamaraca","Itapissuma","Jaboatao dos Guararapes","Moreno","Olinda","Paulista",
            "Pombos","Recife","Sao Lourenco da Mata","Vitoria de Santo Antao"),
  "II"  = c("Bom Jardim","Buenos Aires","Carpina","Casinhas","Cumaru","Feira Nova","Joao Alfredo",
            "Lagoa de Itaenga","Lagoa do Carro","Limoeiro","Machados","Nazare da Mata","Orobo",
            "Passira","Paudalho","Salgadinho","Surubim","Tracunhaem","Vertente do Lerio","Vicencia"),
  "III" = c("Agua Preta","Amaraji","Barreiros","Belem de Maria","Catende","Cortes","Escada",
            "Gameleira","Jaqueira","Joaquim Nabuco","Lagoa dos Gatos","Maraial","Palmares",
            "Primavera","Quipapa","Ribeirao","Rio Formoso","Sao Benedito do Sul",
            "Sao Jose da Coroa Grande","Sirinhaem","Tamandare","Xexeu"),
  "IV"  = c("Agrestina","Alagoinha","Altinho","Barra de Guabiraba","Belo Jardim","Bezerros","Bonito",
            "Brejo da Madre de Deus","Cachoeirinha","Camocim de Sao Felix","Caruaru","Cupira",
            "Frei Miguelinho","Gravata","Ibirajuba","Jatauba","Jurema","Panelas","Pesqueira","Pocao",
            "Riacho das Almas","Saire","Sanharo","Santa Cruz do Capibaribe","Santa Maria do Cambuca",
            "Sao Bento do Una","Sao Caetano","Sao Joaquim do Monte","Tacaimbo",
            "Taquaritinga do Norte","Toritama","Vertentes"),
  "V"   = c("Aguas Belas","Angelim","Bom Conselho","Brejao","Caetes","Calcado","Canhotinho",
            "Capoeiras","Correntes","Garanhuns","Iati","Itaiba","Jucati","Jupi","Lagoa do Ouro",
            "Lajedo","Palmeirina","Paranatama","Saloa","Sao Joao","Terezinha"),
  "VI"  = c("Arcoverde","Buique","Custodia","Ibimirim","Inaja","Jatoba","Manari","Pedra",
            "Petrolandia","Sertania","Tacaratu","Tupanatinga","Venturosa"),
  "VII" = c("Belem de Sao Francisco","Cedro","Mirandiba","Salgueiro","Serrita","Terra Nova","Verdejante"),
  "VIII"= c("Afranio","Cabrobo","Dormentes","Lagoa Grande","Oroco","Petrolina","Santa Maria da Boa Vista"),
  "IX"  = c("Araripina","Bodoco","Exu","Granito","Ipubi","Moreilandia","Ouricuri","Parnamirim",
            "Santa Cruz","Santa Filomena","Trindade"),
  "X"   = c("Afogados da Ingazeira","Brejinho","Carnaiba","Iguaraci","Ingazeira","Itapetim","Quixaba",
            "Santa Terezinha","Sao Jose do Egito","Solidao","Tabira","Tuparetama"),
  "XI"  = c("Betania","Calumbi","Carnaubeira da Penha","Flores","Floresta","Itacuruba",
            "Santa Cruz da Baixa Verde","Sao Jose do Belmonte","Serra Talhada","Triunfo"),
  "XII" = c("Goiana","Alianca","Camutanga","Condado","Ferreiros","Itambe","Itaquitinga",
            "Macaparana","Sao Vicente Ferrer","Timbauba")
)

geres_df <- bind_rows(lapply(names(geres_list), function(g) tibble(geres = g, geres_muni_name = geres_list[[g]])))
stopifnot(nrow(geres_df) == 185, !any(duplicated(geres_df$geres_muni_name)))

macro_map <- c("I" = "Metropolitana", "II" = "Metropolitana", "III" = "Metropolitana", "XII" = "Metropolitana",
               "IV" = "Agreste", "V" = "Agreste",
               "VI" = "Sertao", "X" = "Sertao", "XI" = "Sertao",
               "VII" = "Vale do Sao Francisco e Araripe", "VIII" = "Vale do Sao Francisco e Araripe", "IX" = "Vale do Sao Francisco e Araripe")
geres_df$official_macroregion <- macro_map[geres_df$geres]
stopifnot(!any(is.na(geres_df$official_macroregion)))

# ---- Match against the project's authoritative IBGE municipality lookup ----
lookup <- readRDS(file.path(root, "01_Data/ibge_muni_name_lookup.rds")) |> filter(uf == "PE")
stopifnot(nrow(lookup) == 185)

normalize_name <- function(x) {
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- str_to_lower(x)
  str_replace_all(x, "[^a-z]", "")
}
lookup <- lookup |> mutate(name_norm = normalize_name(name_muni))
geres_df <- geres_df |> mutate(name_norm = normalize_name(geres_muni_name))

# 5 spelling/preposition variants between the GERES source list and the
# project's IBGE lookup (same municipality in every case; verified by IBGE
# code and Google Maps location, not a substantive reclassification):
#   "Gloria de Goita" (GERES doc) = "Glória do Goitá" (IBGE, muni6 260610)
#   "Itamaraca"        (GERES doc) = "Ilha de Itamaracá" (IBGE, muni6 260760)
#   "Sao Caetano"      (GERES doc) = "São Caitano" (IBGE, muni6 261310)
#   "Belem de Sao Francisco" (GERES doc) = "Belém do São Francisco" (IBGE, muni6 260160)
#   "Iguaraci"         (GERES doc) = "Iguaracy" (IBGE, muni6 260690)
spelling_fix <- c(gloriadegoita = "gloriadogoita", itamaraca = "ilhadeitamaraca",
                   saocaetano = "saocaitano", belemdesaofrancisco = "belemdosaofrancisco",
                   iguaraci = "iguaracy")
geres_df <- geres_df |> mutate(name_norm = if_else(name_norm %in% names(spelling_fix), spelling_fix[name_norm], name_norm))

matched <- geres_df |> left_join(lookup |> select(muni6, name_muni, name_norm), by = "name_norm")
if (any(is.na(matched$muni6))) {
  message("Unmatched:"); print(as.data.frame(matched |> filter(is.na(muni6))))
  stop("STOP: unresolved municipality name mismatch between the official GERES list and the project's IBGE lookup -- investigate before building the stratum lookup.")
}
if (n_distinct(matched$muni6) != 185) stop("STOP: matched municipality set is not exactly the 185 PE administrative units.")
message("CONFIRMED: all 185 official GERES municipality entries matched 1:1 to the project's IBGE muni6 lookup.")

# ---- Final 5-way model stratum ----------------------------------------
RECIFE_MUNI6 <- matched$muni6[matched$name_norm == normalize_name("Recife")]
stopifnot(length(RECIFE_MUNI6) == 1)

matched <- matched |>
  mutate(final_model_stratum = case_when(
    muni6 == RECIFE_MUNI6 ~ "1_Recife",
    official_macroregion == "Metropolitana" ~ "2_Metropolitana_remainder",
    official_macroregion == "Agreste" ~ "3_Agreste",
    official_macroregion == "Sertao" ~ "4_Sertao",
    official_macroregion == "Vale do Sao Francisco e Araripe" ~ "5_Vale_Sao_Francisco_Araripe",
    TRUE ~ NA_character_
  ))
stopifnot(!any(is.na(matched$final_model_stratum)))

# ---- Population: most recent (last observed week) population per muni6 ----
# from the already-built, validated PE municipality-week panel (Section 1 of
# the spatial-turnover diagnostic). Fernando de Noronha is a federal state
# district, not a municipality, and is absent from the case panel and from
# the national demography series used elsewhere in this project -- flagged
# with NA population, not imputed.
pe_panel <- readRDS(file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic/pe_municipality_week_panel.rds"))
pop_latest <- pe_panel |> filter(week_start == max(week_start)) |> select(muni6, population)

# ---- Section 2: four-macroregion sensitivity partition (prepared, NOT fit) ----
# Same official macroregions, WITHOUT separating Recife. Saved as an extra
# column for later use only, per instruction not to fit this yet.
four_stratum_map <- c("Metropolitana" = "A_Metropolitana", "Agreste" = "B_Agreste",
                       "Sertao" = "C_Sertao", "Vale do Sao Francisco e Araripe" = "D_Vale_Sao_Francisco_Araripe")
matched$four_stratum_sensitivity <- four_stratum_map[matched$official_macroregion]
stopifnot(!any(is.na(matched$four_stratum_sensitivity)))

stratum_lookup <- matched |>
  left_join(pop_latest, by = "muni6") |>
  transmute(
    muni6, municipality_name = name_muni, geres, official_macroregion,
    final_model_stratum, four_stratum_sensitivity, population
  ) |>
  arrange(final_model_stratum, desc(population))

write_csv(stratum_lookup, file.path(table_dir, "PE_spatial_stratum_lookup.csv"))
message("\n[saved] ", file.path(table_dir, "PE_spatial_stratum_lookup.csv"))

message("\n=== Stratum summary ===")
print(as.data.frame(stratum_lookup |> group_by(final_model_stratum) |>
                       summarise(n_municipalities = n(),
                                 n_missing_population = sum(is.na(population)),
                                 total_population_2025 = sum(population, na.rm = TRUE), .groups = "drop")))

message("\n=== Municipalities with missing population (expected: Fernando de Noronha only) ===")
print(as.data.frame(stratum_lookup |> filter(is.na(population)) |> select(muni6, municipality_name, geres, final_model_stratum)))
