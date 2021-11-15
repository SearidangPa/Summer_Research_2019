import numpy as np 
import random 
import cmath 
import math
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
pi = cmath.pi

def main(N, r, l, equal_dist, Folder_name):    
    file1 = open("Log_"+ str(Folder_name)+".txt",'w') 
    mb = math.ceil(N/r) 
    num_iterations = 2
    for x in range(num_iterations):
        print("performing experiment number: " + str(x+1))

        ListK = [0, *generate_random_insertion_or_deletion(mb, l, equal_dist), mb]
        [Shift_at_each_interval, Insert_or_Delete] = Get_Shift_at_each_interval(l)
        
        Y = np.arange(0, N, 1)
        f = probability_distr_with_insertion_deletion(Y, N, mb, r, ListK, Shift_at_each_interval)
        plot_figure(Y, f, x, Folder_name, ListK, Insert_or_Delete)
        write_log(file1, ListK, Shift_at_each_interval, x) 

    file1.close()
       
def plot_figure(Y, f, x, Folder_name, ListK, Insert_or_Delete):
    for y in Y:
        if  (y*r) % N <= math.floor(r/2) or (y*r) % N >= N-math.floor(r/2):
            prob_new_ys = new_good_y(f, y, l)
            plt.plot(Y[y], prob_new_ys, 'ko')

            plt.plot(Y[y], f[y], 'bo')
            plt.plot(Y[y-1], f[y-1], 'go')
            plt.plot(Y[y+1], f[y+1], 'yo')
            if l >=2:
                plt.plot(Y[y-2], f[y-2], 'co')
                plt.plot(Y[y+2], f[y+2], 'mo')
            if l >=3:
                plt.plot(Y[y-3], f[y-3], color = 'pink', marker = 'o')
                plt.plot(Y[y+3], f[y+3], color = 'indigo', marker ='o')   
            if l >=4:
                plt.plot(Y[y-4], f[y-4], color = 'peru', marker = 'o')
                plt.plot(Y[y+4], f[y+4], color = 'sandybrown', marker ='o')
            if l >=5:
                plt.plot(Y[y-5], f[y-5], color = 'darkblue', marker = 'o')
                plt.plot(Y[y+5], f[y+5], color = 'darkslategray', marker ='o')
            if l >=6:
                plt.plot(Y[y-6], f[y-6], color = 'orange', marker = 'o')
                plt.plot(Y[y+6], f[y+6], color = 'olive', marker ='o') 
            if l >=7:
                plt.plot(Y[y-7], f[y-7], color = 'mediumspringgreen', marker = 'o')
                plt.plot(Y[y+7], f[y+7], color = 'gold', marker ='o') 
    
    plot_insertion_deletion(ListK, Insert_or_Delete)
    plt.savefig(Folder_name + "/experiment_number" + str(x+1) + ".png") 
    plt.clf()   #clear figure 


def plot_insertion_deletion(ListK, Insert_or_Delete):
    for i in range(1, len(ListK)-1):
        if Insert_or_Delete[i] > 0:
            plt.plot(ListK[i]*r, -1/(100*r), 'r*')
        else:
            plt.plot(ListK[i]*r, -1/(100*r), 'ro')


def write_log(file1, ListK, Shift_at_each_interval, x):
    #----------------write to file to keep a log of the experiments------------------------ 
    file1.write('experiment number: ' + str(x+1))
    file1.write(': N: ' + str(N) + ', r: '+ str(r) + ', l: '+ str(l))
    file1.write('\nListK: ') 
    for i in ListK:
        file1.write(str(i)+' ') 

    file1.write('\nShift: ')
    for i in Shift_at_each_interval:
        file1.write(str(math.floor(i))+' ')
    file1.write("\n\n")



#----return a list of numK elements sampled from a range (1, mb). Each acts as a 
#boundary for each interval of insertions or deletions (just one insertion at a time for now)
def generate_random_insertion_or_deletion(mb, l, equal_dist):
    #start_sample = math.floor(3*mb/4)
    #end_sample = math.floor(3*mb/4)+l
    start_sample = 1
    end_sample = mb
    ListK = sorted(random.sample(range(start_sample, end_sample), l)) 
    if equal_dist == True:
        for i in range(0, l):
            ListK[i] = math.ceil((i+1)*(mb/(l+1)))
    return ListK


def Get_Shift_at_each_interval(l):
    Insert_or_Delete = (l+1)*[0]
    if random.randint(0, 9) < 5:
        for i in range(1, l+1):
            Insert_or_Delete[i] = 1
    else:
        for i in range(1, l+1):
            Insert_or_Delete[i] = -1
    Shift_at_each_interval = np.cumsum(Insert_or_Delete)
    return [Shift_at_each_interval, Insert_or_Delete]


def new_good_y (f, y, l):
    result = f[y]
    for i in range(1, l+1):
        result = result + f[y-i] + f[y+i]
    return result 

#----------The probability distribution when there are insertions or deletions-------------
def probability_distr_with_insertion_deletion(Y, N, mb, r, ListK, Shift_at_each_interval):
    f = np.zeros(N)         #initialize the result of the probability for each y 
    for y in Y:             #for each y, calculate the probability 
        z = complex(0,0)    #store the total sum of the complex number in z 
        for i in range(len(ListK)-1):
            temp_sum = complex(0,0)     #store a temporary sum between each interval 
            for m in range(ListK[i], ListK[i+1]):
                temp_sum = temp_sum + cmath.exp(2*pi*1j*m*r*y/N)    
            z = z + temp_sum * cmath.exp(2*pi*y*1j*Shift_at_each_interval[i]/N)
        f[y] = abs(z)**2/(mb*N)
    return f


if __name__ == "__main__":
    equal_dist = False      
    Folder_name = "Plot1"                            
    
    N = 4000
    r = 100
    l = 3
    main(N, r, l,  equal_dist, Folder_name)







