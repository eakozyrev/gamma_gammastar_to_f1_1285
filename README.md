# gamma_gammastar_to_f1_1285
scripts to analyse the process of two gamma fusion into f1_1285 meson

# Обработка root файлов
root -l make.C

# Аппроксимация спекторов 2пи эта (итоговый файл ../eef1/results/nevents.dat)
root -l
.L alpha.C+
fitall()

# Вычисление сечения
root -l
.L alpha1.C+
cross_section()


# Разделение m = 1 и m = 0 
root -l
.L alpha1.C+
script_minim()