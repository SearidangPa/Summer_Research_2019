import numpy as np 
import random 
import cmath 
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors

pi = cmath.pi

def main(N, r, l, file1): 
    mb = math.ceil(N/r)    
    f1 = open(file1+ '_log.txt', 'w')  
    num_iterations = 50
   
    for x in range(num_iterations):
        print("performing experiment number: " + str(x+1))
        start_sample = 40
        sample = generate_random_insertion_or_deletion(start_sample+1, mb-2, l)
        Listk1 = [0, *sample, mb]
        Listk2 = [0, start_sample, *sample, mb]

        [Shift1, Ins_Del1] = Get_Shift1(l)
        if Ins_Del1[1] ==-1:
            Ins_Del2 = [*Ins_Del1, -1]
        else:
            Ins_Del2 = [*Ins_Del1, 1]
        Shift2 = np.cumsum(Ins_Del2)

        print(Listk1)
        print(Shift1)
        print(Listk2)
        print(Shift2)

        y = np.arange(0, N, 1)
        function1 = probability_distr_with_insertion_deletion(y, N, mb, r, Listk1, Shift1)
        function2 = probability_distr_with_insertion_deletion(y, N, mb, r, Listk2, Shift2)

        prob_y_worst1 = 1
        prob_y_worst2 = 1
        for t in range(1, N):
            if  (t*r) % N <= math.floor(r/2) or (t*r) % N >= N-math.floor(r/2):                                                  
                
                prob_new_y1 = new_good_y(function1, t, l)
                plt.plot(y[t], prob_new_y1, 'ko')

                prob_new_y2 = new_good_y(function2, t, l+1)
                plt.plot(y[t], prob_new_y2, 'mo')
                
                plt.plot(y[t], function1[t], 'bo')
                plt.plot(y[t], function2[t], 'go')

                plt.plot(y[t], function1[t-1], 'yo')
                plt.plot(y[t], function2[t-1], 'co')

                if prob_new_y1 < prob_y_worst1:
                    prob_y_worst1 = prob_new_y1
                if prob_new_y2 < prob_y_worst2:
                    prob_y_worst2 = prob_new_y2
                
        for i in range(1, len(Listk2)-1):
            if Ins_Del2[i] > 0:
                plt.plot(y[Listk2[i]]*r, -1/(100*r), 'r*')
            else:
                plt.plot(y[Listk2[i]]*r, -1/(100*r), 'ro')
                
        if prob_y_worst1 > prob_y_worst2:
            diff = round((prob_y_worst1 - prob_y_worst2), 8)
            plt.xlabel("prob_y_worst_1 - prob_y_worst_2 = "+str(diff))
        
        plt.savefig(file1 + "/experiment_number" + str(x+1) + ".png") 
        plt.clf()   #clear figure 

        #----------------write to file to keep a log of the experiments------------------------ 
        f1.write('experiment number: ' + str(x+1))
        f1.write(': N: ' + str(N) + ', r: '+ str(r) + ', l: '+ str(l))
        f1.write('\nListk1: ') 
        for i in Listk1:
            f1.write(str(i)+' ') 

        f1.write('\nShift: ')
        for i in Shift1:
            f1.write(str(math.floor(i))+' ')

        f1.write("\n\n")
        f1.flush()
    f1.close()

       
def Get_Shift1(l):
    Ins_Del1 = (l+1)*[0]
    if random.randint(0, 9) < 5:
        for i in range(1, l+1):
            Ins_Del1[i] = 1
    else:
        for i in range(1, l+1):
            Ins_Del1[i] = -1

    Shift1 = np.cumsum(Ins_Del1)
    
    return [Shift1, Ins_Del1]

#----return a list of numK elements sampled from a range (1, mb). Each acts as a 
#boundary for each interval of insertions or deletions (just one insertion at a time for now)
def generate_random_insertion_or_deletion(start_sample, mb, l):
    Listk1 = sorted(random.sample(range(start_sample, mb), l)) 
    return Listk1
   
def new_good_y (f, t, l):
    result = f[t]
    for i in range(1, l+1):
        result = result + f[t-i] + f[t+i]
    return result 

#----------The probability distribution when there are insertions or deletions-------------
def probability_distr_with_insertion_deletion(y, N, mb, r, Listk1, Shift1):
    f = np.zeros(N)         #initialize the result of the probability for each y 
    for t in y:             #for each y, calculate the probability 
        z = complex(0,0)    #store the total sum of the complex number in z 
        for i in range(len(Listk1)-1):
            temp_sum = complex(0,0)     #store a temporary sum between each interval 
            for m in range(Listk1[i], Listk1[i+1]):
                temp_sum = temp_sum + cmath.exp(2*pi*1j*m*r*t/N)    
            z = z + temp_sum * cmath.exp(2*pi*t*1j*Shift1[i]/N)
        f[t] = abs(z)**2/(mb*N)
    return f


if __name__ == "__main__":
    #folder to save into                                 
    file1 = 'Test'

    N = 10000
    r = 50
    l = 2

    print('N = ' + str(N) + ', r= '+ str(r) + ', l= '+ str(l))
    main(N, r, l, file1)


    







