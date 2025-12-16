################################################################################
# MLN-SEM Shiny App - 파일 확인 및 테스트 스크립트
#
# 앱 실행 전에 모든 필수 파일이 있는지 확인합니다.
################################################################################

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  MLN-SEM SHINY APP - 파일 확인 및 테스트                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# 필수 파일 목록
required_files <- list(
  "핵심 앱 파일" = c(
    "app.R"
  ),
  "분석 프레임워크" = c(
    "MLN_SEM_v2_1_FINAL.R",
    "MLN_SEM_Visualization_FIXED.R",
    "MLN_SEM_Extended_Network_Visualization_v2.R",
    "MLN_SEM_Module_Assignments_Simple.R"
  ),
  "예시 데이터" = c(
    "example_clinical.csv",
    "example_metabolome.csv",
    "example_microbiome.csv"
  )
)

# 선택 파일 목록
optional_files <- list(
  "설치/배포 스크립트" = c(
    "install_packages.R",
    "deploy_shinyapps.R"
  ),
  "문서" = c(
    "README.md",
    "한글_사용가이드.md",
    "QUICKSTART.md",
    "DEPLOYMENT_GUIDE.md"
  )
)

# 파일 확인 함수
check_files <- function(file_list, category_name, required = TRUE) {
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("%s: %s\n", category_name, if(required) "필수" else "선택"))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"))
  
  all_present <- TRUE
  
  for (file in file_list) {
    exists <- file.exists(file)
    
    if (exists) {
      size <- file.info(file)$size
      size_kb <- round(size / 1024, 1)
      cat(sprintf("  ✓ %s (%.1f KB)\n", file, size_kb))
    } else {
      cat(sprintf("  ✗ %s [누락]\n", file))
      all_present <- FALSE
    }
  }
  
  cat("\n")
  return(all_present)
}

# 작업 디렉토리 확인
cat("현재 작업 디렉토리:\n")
cat(sprintf("  %s\n\n", getwd()))

# 필수 파일 확인
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  필수 파일 확인                                                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

all_required_present <- TRUE

for (category in names(required_files)) {
  present <- check_files(required_files[[category]], category, required = TRUE)
  all_required_present <- all_required_present && present
}

# 선택 파일 확인
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  선택 파일 확인                                                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

for (category in names(optional_files)) {
  check_files(optional_files[[category]], category, required = FALSE)
}

# 최종 결과
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║  확인 결과                                                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

if (all_required_present) {
  cat("  ✓✓✓ 모든 필수 파일이 있습니다! ✓✓✓\n\n")
  cat("앱을 실행할 준비가 되었습니다:\n\n")
  cat("  library(shiny)\n")
  cat("  runApp('app.R')\n\n")
  
  # 예시 데이터 로드 테스트
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("예시 데이터 로드 테스트 중...\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  
  tryCatch({
    clinical <- read.csv("example_clinical.csv", stringsAsFactors = FALSE)
    metabolome <- read.csv("example_metabolome.csv", stringsAsFactors = FALSE)
    microbiome <- read.csv("example_microbiome.csv", stringsAsFactors = FALSE)
    
    cat(sprintf("  ✓ Clinical data:   %d × %d\n", nrow(clinical), ncol(clinical)))
    cat(sprintf("  ✓ Metabolome data: %d × %d\n", nrow(metabolome), ncol(metabolome)))
    cat(sprintf("  ✓ Microbiome data: %d × %d\n", nrow(microbiome), ncol(microbiome)))
    cat("\n  → 예시 데이터가 정상적으로 로드됩니다!\n\n")
    
  }, error = function(e) {
    cat(sprintf("  ✗ 에러: %s\n\n", e$message))
  })
  
  return(invisible(TRUE))
  
} else {
  cat("  ✗✗✗ 일부 필수 파일이 누락되었습니다! ✗✗✗\n\n")
  cat("다음 작업을 수행하세요:\n\n")
  cat("1. 누락된 파일들을 다운로드\n")
  cat("2. 모든 파일을 같은 디렉토리에 배치\n")
  cat("3. 이 스크립트를 다시 실행\n\n")
  cat("도움말:\n")
  cat("  - README.md 또는 FILE_STRUCTURE.md 참고\n")
  cat("  - QUICKSTART.md에서 빠른 시작 가이드 확인\n\n")
  
  return(invisible(FALSE))
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("필요한 패키지를 설치하려면:\n")
cat("  source('install_packages.R')\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
