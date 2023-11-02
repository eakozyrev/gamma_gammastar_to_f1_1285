import os

max_ = 25
min_ = 0.02
Q2_ = 0.02

file_ = open("cross_secs.txt","w")

while Q2_ <= max_:
    cmd_ = './test {}'.format(-Q2_)
    os.system(cmd_)
    file0 = open("cross.txt","r")
    cr = file0.read()
    file_.write('{:.2f} {:.8f} \n'.format(float(Q2_),float(cr)))
    Q2_ = Q2_ + min_
    file0.close()

    
