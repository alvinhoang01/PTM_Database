del dist\*.whl
del /s /q build\*

python setup.py bdist_wheel

pip uninstall qcmspycloud --yes

pip install dist\qcmspycloud-0.2.0-py3-none-any.whl

REM qcmspycloud "F:\OneDrive - Johns Hopkins\project\Group\Liwei\core_fuc\output\Fetuin_reproducibility"