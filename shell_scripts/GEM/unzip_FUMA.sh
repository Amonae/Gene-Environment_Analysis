srun -p sunlab --pty bash

hormone="SHBG"
interaction=(BMI T2D CAD)
sex=(Female Male)

for x in "${interaction[@]}";
do for y in "${sex[@]}";
do dir="/projects/sunlab/Students_Work/Amonae_work/GEM_$hormone/$x/$y/FUMA_results";
cd $dir;
unzip \*.zip
rm *.zip
done;
done

