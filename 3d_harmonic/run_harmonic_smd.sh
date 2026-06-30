i=$1
strain=$2
temp=$3
python3 generate_smd_harmonic.py 50 --strain $strain
python3 smd_harmonic.py \
    --temperature $temp \
    --ensemble nvt \
    --dimensionality 3 \
    --total_integrated_atoms 50 \
    --log_name "3d_smd_rc/log.T${temp}_s${strain}_${i}" \
    --strain $strain \
    --data_file final.lmp && \
mpiexec -np 1 /usr/workspace/potter21/lmp_lassen -in in.harmonic_smd -screen out
