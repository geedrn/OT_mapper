@echo off
REM =============================================================================
REM OT Mapper Shiny App Launcher (Windows)
REM =============================================================================
REM This script helps users launch the Shiny app easily on Windows.
REM
REM Usage: Double-click run_app.bat or run from command prompt
REM =============================================================================

echo ==========================================
echo OT Mapper Shiny App Launcher
echo ==========================================
echo.

REM Change to script directory
cd /d "%~dp0"

REM Check if R is installed
where Rscript >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] R is not installed or not in PATH
    echo.
    echo Please install R first:
    echo   Download from https://cran.r-project.org/
    echo   Make sure to add R to PATH during installation
    echo.
    pause
    exit /b 1
)

echo [OK] R is installed
echo.

REM Check if bedtools is installed
where bedtools >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [WARNING] bedtools is not installed
    echo.
    echo bedtools is required for gene annotation. The app will work but
    echo annotation features will be disabled.
    echo.
    echo To install bedtools, download from:
    echo   https://bedtools.readthedocs.io/en/latest/content/installation.html
    echo.
    set /p continue="Continue anyway? (y/n): "
    if /i not "%continue%"=="y" exit /b 1
) else (
    echo [OK] bedtools is installed
)
echo.

REM Check for required R packages
echo Checking required R packages...
Rscript -e "required_packages <- c('shiny', 'shinydashboard', 'DT', 'httr', 'readr', 'dplyr', 'stringr', 'GenomicRanges'); missing_packages <- c(); for (pkg in required_packages) { if (!requireNamespace(pkg, quietly = TRUE)) { missing_packages <- c(missing_packages, pkg) } }; if (length(missing_packages) > 0) { cat('Installing missing packages...\n'); cran_packages <- missing_packages[missing_packages != 'GenomicRanges']; if (length(cran_packages) > 0) { install.packages(cran_packages, repos='https://cran.rstudio.com/', quiet=TRUE) }; if ('GenomicRanges' %in% missing_packages) { if (!requireNamespace('BiocManager', quietly = TRUE)) { install.packages('BiocManager', repos='https://cran.rstudio.com/', quiet=TRUE) }; BiocManager::install('GenomicRanges', quiet=TRUE) }; cat('Packages installed successfully!\n') } else { cat('All required packages are installed.\n') }"
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Failed to check/install R packages
    pause
    exit /b 1
)

echo.
echo [OK] All dependencies are ready!
echo.

REM Check if annotation data files exist
if not exist "..\R\scripts\data\UCSC_exons_modif_canonical.bed" (
    echo [WARNING] Annotation database files not found
    echo.
    echo The app will work, but gene annotation features will be disabled.
    echo Annotation files should be in: R\scripts\data\
    echo.
    set /p continue="Continue anyway? (y/n): "
    if /i not "%continue%"=="y" exit /b 1
) else (
    echo [OK] Annotation database files found
)

echo.
echo ==========================================
echo Starting Shiny App...
echo ==========================================
echo.
echo The app will open in your default web browser.
echo If it doesn't open automatically, go to: http://localhost:3838
echo.
echo Press Ctrl+C to stop the app
echo.

REM Run the Shiny app
Rscript -e "shiny::runApp(port=3838, launch.browser=TRUE)"

pause

