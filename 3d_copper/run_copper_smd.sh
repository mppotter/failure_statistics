i=$1
strain=$2
temp=$3
python3 generate_smd_Cu.py 50 --strain $strain
python3 smd_copper.py \
    --temperature $temp \
    --ensemble nvt \
    --dimensionality 3 \
    --total_integrated_atoms 50 \
    --log_name "copper_smd/log.T${temp}_s${strain}_${i}" \
    --strain $strain \
    --meam_file Cu.meam \
    --library_file library.meam \
    --table_file harmonic_potential.txt \
    --break_r_file Cu_breaking_rs.txt \
    --data_file final.lmp && \
flux run -n 1 /usr/workspace/potter21/lammps/build/lmp -in in.copper_smd
