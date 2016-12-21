@echo off
for %%a in ( *.job ) do (
    java -jar ..\dist\SRA_DNLD_MGR.jar %%a
    if ERRORLEVEL 1 (
        echo error job %%a
        exit /b 1
    )
)
