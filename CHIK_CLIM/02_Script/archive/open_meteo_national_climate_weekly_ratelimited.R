# ===========================================================================
# [ARCHIVED] open_meteo_national_climate_weekly_ratelimited.R
# (원래 이름: 04_fetch_national_climate_weekly.R)
#
# 왜 archive로 옮겼나
# ------------------
# muni 5,122개 x 10년치를 Open-Meteo Archive API로 받으려다 HTTP 429로
# 막힘. 문서상 무료 한도(600 call/min, 10,000 call/day)는 넉넉해 보였지만,
# 지점 수 x 연도 수 x 변수 수를 내부적으로 call로 환산하는 것으로 보여
# batch를 아무리 늦춰도 뚫리지 않았음. 무료 API key도 존재하지 않음(무료
# tier = key 없는 anonymous tier 그 자체).
#
# 대체
# ----
# 30_climate_covariates_dlnm/04_fetch_era5land_temperature.R +
# 05_fetch_era5land_precip.R + 06_build_national_climate_weekly.R 로
# 대체됨 — Copernicus CDS에서 브라질 전역 grid를 통짜로 받아(지점별 요청이
# 아니라 국가 단위 요청이라 rate-limit 문제 자체가 없음) muni centroid에서
# 로컬로 추출하는 방식. Open-Meteo Archive API도 결국 ERA5-Land를 감싼
# 것이라 데이터 자체는 동일.
#
# 목적 (원래 스크립트 설명, 참고용)
# ----

# chik_brazil_muni_week_2015_2024.rds 에 있는 모든 municipality(muni6)에
# 대해, Open-Meteo Archive API에서 muni centroid 기준 daily 기후를 받아서
# case panel과 동일한 epi-week(week_start, 일요일 시작) 단위로 집계한다.
#
# 규모
# ----
# ~5,152개 municipality. Open-Meteo Archive API는 위경도를 comma로 여러 개
# 넘기면 한 번의 호출로 여러 지점을 반환하는 batch 조회를 지원한다(실측 확인
# 완료). batch_size=25 기준 약 200회 호출, 전체 약 45~70분 소요 예상.
#
# 재시작 가능성 (checkpointing)
# ------------------------------
# 외부 API를 200회 가까이 호출하는 장시간 작업이므로, batch 하나가 끝날 때마다
# 결과를 01_Data/climate_cache/batch_####.rds 로 저장한다. 스크립트를 다시
# 실행하면 이미 캐시된 batch는 건너뛰고 이어서 진행한다 (처음부터 다시 받지
# 않음).
#
# timezone
# --------
# Open-Meteo의 timezone=auto 대신 전 지점에 timezone="America/Sao_Paulo"를
# 고정한다. auto를 쓰면 지역별 시간대(예: 서부 일부 UTC-4/-5)에 따라 주
# 경계(week_start)가 지점마다 미묘하게 어긋날 수 있기 때문.
#
# 출력
# ----
#   01_Data/chik_muni_week_climate_2015_2024.rds
#     columns: muni6, week_start, PRCP (주간 합), Tmin/Tmax/Tmean (주간
#              최저/최고/평균), n_days (집계에 쓰인 일수, 정상은 7)
# ===========================================================================

for (p in c("here", "httr2", "dplyr", "lubridate", "tibble", "readr", "purrr")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(here); library(httr2); library(dplyr); library(lubridate)
  library(tibble); library(readr); library(purrr)
})

fetch_climate_batch <- function(munis_batch, start_date, end_date, timezone) {
  lat_str <- paste(munis_batch$lat, collapse = ",")
  lon_str <- paste(munis_batch$lon, collapse = ",")

  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    "latitude=", lat_str, "&longitude=", lon_str,
    "&start_date=", start_date, "&end_date=", end_date,
    "&daily=precipitation_sum,temperature_2m_min,temperature_2m_max,temperature_2m_mean",
    "&timezone=", utils::URLencode(timezone, reserved = TRUE)
  )

  resp <- httr2::request(url) |>
    httr2::req_timeout(180) |>
    httr2::req_retry(max_tries = 6, backoff = \(i) min(20 * i, 90)) |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  if (length(resp) != nrow(munis_batch)) {
    stop(sprintf(
      "Open-Meteo returned %d locations but batch had %d municipalities.",
      length(resp), nrow(munis_batch)
    ))
  }

  purrr::map2_dfr(resp, munis_batch$muni6, function(loc, muni6) {
    tibble::tibble(
      muni6     = muni6,
      date      = as.Date(unlist(loc$daily$time)),
      PRCP_day  = as.numeric(unlist(loc$daily$precipitation_sum)),
      Tmin_day  = as.numeric(unlist(loc$daily$temperature_2m_min)),
      Tmax_day  = as.numeric(unlist(loc$daily$temperature_2m_max)),
      Tmean_day = as.numeric(unlist(loc$daily$temperature_2m_mean))
    )
  })
}

fetch_national_climate_weekly <- function(
    muni6_vec,
    week_start_vec,
    centroids,
    batch_size = 50,
    cache_dir  = here::here("01_Data/climate_cache"),
    timezone   = "America/Sao_Paulo",
    pause_sec  = 20
) {
  week_start_vec <- sort(unique(as.Date(week_start_vec)))
  start_date <- min(week_start_vec)
  end_date   <- max(week_start_vec) + 6L

  munis <- centroids |>
    dplyr::filter(muni6 %in% muni6_vec) |>
    dplyr::arrange(muni6)

  unmatched <- setdiff(muni6_vec, munis$muni6)
  if (length(unmatched) > 0) {
    message(sprintf(
      "[warn] %d muni6 in the case panel have no centroid and will be skipped: %s%s",
      length(unmatched),
      paste(head(unmatched, 10), collapse = ", "),
      if (length(unmatched) > 10) ", ..." else ""
    ))
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  batch_id <- ceiling(seq_len(nrow(munis)) / batch_size)
  batches  <- split(munis, batch_id)
  n_batch  <- length(batches)

  message(sprintf(
    "[fetch] %d municipalities, %d batches of <=%d, %s to %s",
    nrow(munis), n_batch, batch_size, start_date, end_date
  ))

  for (i in seq_along(batches)) {
    cache_file <- file.path(cache_dir, sprintf("batch_%04d.rds", i))
    if (file.exists(cache_file)) next

    b <- batches[[i]]
    daily <- tryCatch(
      fetch_climate_batch(b, start_date, end_date, timezone),
      error = function(e) {
        message(sprintf("[batch %d/%d] FAILED: %s", i, n_batch, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(daily)) next

    saveRDS(daily, cache_file)
    message(sprintf("[batch %d/%d] done (%d municipalities)", i, n_batch, nrow(b)))
    Sys.sleep(pause_sec)
  }

  cache_files <- list.files(cache_dir, pattern = "^batch_\\d{4}\\.rds$", full.names = TRUE)
  fetched_batches <- length(cache_files)
  if (fetched_batches < n_batch) {
    stop(sprintf(
      "Only %d/%d batches cached — some failed after retries. Re-run this script to resume.",
      fetched_batches, n_batch
    ))
  }

  message("[aggregate] combining daily data from ", fetched_batches, " cached batches ...")
  daily_all <- purrr::map_dfr(cache_files, readRDS)

  daily_all <- daily_all |>
    dplyr::mutate(week_start = lubridate::floor_date(date, unit = "week", week_start = 7))

  weekly <- daily_all |>
    dplyr::filter(week_start %in% week_start_vec) |>
    dplyr::group_by(muni6, week_start) |>
    dplyr::summarise(
      PRCP   = sum(PRCP_day, na.rm = TRUE),
      Tmin   = min(Tmin_day, na.rm = TRUE),
      Tmax   = max(Tmax_day, na.rm = TRUE),
      Tmean  = mean(Tmean_day, na.rm = TRUE),
      n_days = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(muni6, week_start)

  weekly
}

if (sys.nframe() == 0) {
  # ============================================================
  # CONFIG — set N_MUNI_LIMIT to a small number (e.g. 5) for a
  # quick smoke test before committing to the full ~1hr national run.
  # ============================================================
  N_MUNI_LIMIT <- NULL
  BATCH_SIZE   <- 50
  # ============================================================

  panel <- readRDS(here::here("01_Data/chik_brazil_muni_week_2015_2024.rds"))
  centroids <- readRDS(here::here("01_Data/ibge_muni_centroids.rds"))

  target_munis <- sort(unique(panel$muni6))
  if (!is.null(N_MUNI_LIMIT)) {
    target_munis <- head(target_munis, N_MUNI_LIMIT)
    message("[smoke test] limited to ", N_MUNI_LIMIT, " municipalities")
  }

  weekly_climate <- fetch_national_climate_weekly(
    muni6_vec      = target_munis,
    week_start_vec = panel$week_start,
    centroids      = centroids,
    batch_size     = BATCH_SIZE
  )

  out_path <- here::here("01_Data/chik_muni_week_climate_2015_2024.rds")
  saveRDS(weekly_climate, out_path)
  readr::write_csv(weekly_climate, sub("\\.rds$", ".csv", out_path))

  message(sprintf(
    "[save] %s (%d rows, %d municipalities, %d weeks)",
    out_path, nrow(weekly_climate),
    dplyr::n_distinct(weekly_climate$muni6),
    dplyr::n_distinct(weekly_climate$week_start)
  ))
}
