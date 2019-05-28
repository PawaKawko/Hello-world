#The main idea here that we divided the space in which we are working in N^(2n+1) dots.
#To each start conditions we use the method of LS to find the minimum of the sum of squares.
#At the end we choose the minimum with the smallest error
#Additionally it will calculate the medium width of phase curve and cut not suitable data N times
#



"""==========================================="""
"""IMPORTING LIBRUARIES"""
"""==========================================="""

import scipy.optimize as spo     #for the method of LS
import numpy as np               #for math stuff
import matplotlib.pyplot as plt  #for plotting
import time                      #to know time of calculations
import tkinter as tnk            #graphic interface
import os                        #to work with directories

"""==========================================="""
"""CONSTANTS (DO NOT CHANGE!!!)"""
"""==========================================="""

N = 3           #dots in each direction (more) 

"""==========================================="""
"""TRIGONOMETRIC POLYNOMIAL FUNCTIONS"""
"""==========================================="""

def sin(t, pp, n):                  #approximation of function by Fourie series (t -> x_data, pp - parameters)
    x = np.zeros(len(t))             #creating array
    for i in range(n):
        x += pp[2*i+1]*np.sin(2*np.pi*t/pp[0]*(i+1)+pp[2*i+2])
    return x

def sin1(t, pp, n, j):              #approximation of function by Fourie series and get y[j]
    y = np.zeros(len(t))
    for i in range(n):
        y += pp[2*i+1]*np.sin(2*np.pi*t/pp[0]*(i+1)+pp[2*i+2])
    return y[j]

"""==========================================="""
"""READING DATA FROM FILE"""
"""==========================================="""

def read_data(name, ftype):
    with open(name + ftype, 'r') as file:   #data
        x, y, y_err = [], [], []             #set arrays
        lines = file.readlines()             #lines - list; read() - not work
        for i in lines: 
            data = i.split()                 #split into words because of spaces
            x.append(float(data[0]))
            y.append(float(data[1]))
            y_err.append(float(data[2]))
        x, y, y_err = np.array(x), np.array(y), np.array(y_err)      #to make arrays more cool and suitable with method of LS
    Number_of_elements = len(x)              #number of observations or data we have
    Number_of_elements_0 = len(x)  
    return x, y, y_err, Number_of_elements, Number_of_elements_0

"""==========================================="""
"""CALCULATING PRESIZE VALUE OF PERIOD"""
"""==========================================="""

def becoming_perfect(bounds, N, x, y, y_err, name, I, Rep):
    T_max, T_min, A_max, A_min, F_max, F_min = bounds[0:6]
    d_T = (T_max - T_min)/N
    d_A = (A_max - A_min)/N
    d_F = (F_max - F_min)/N
    
    p0 = []             #start conditions  
    n = 2
    ans_ideal = np.zeros(2*n+1)              
    Error_0 = 1000      #just for start value
    
    for j_1 in range(2*N):            #we try each start sondition (Period is more important => 2*N)
        for j_2 in range(N):
            for j_3 in range(N):
                for j_4 in range(N):
                    for j_5 in range(N):
                        p0.clear()
                        p0 = [d_T*(j_1 + 1) + T_min, d_A*(j_2 + 1) + A_min, d_F*(j_3 + 1) + F_min, d_A*(j_4 + 1) + A_min, d_F*(j_5 + 1) + F_min]
                                    
                        fun = lambda pp: (y - sin(x, pp, n))/y_err       #core of least squares
                        ans = spo.leastsq(fun, p0, full_output=1)
                        sigma = np.sum((y - sin(x, ans[0], n))**2)/len(x)
                        error = np.sqrt(np.diag(ans[1]*sigma))
                        
                        if np.sum(error) < Error_0:  #if we find better minimum
                            T_ideal = ans[0][0]
                            Error_0 = np.sum(error)
                            error_T = error[0]
                            for i in range(2*n+1):
                                ans_ideal[i]  = ans[0][i]
                               
    order_Error = -int(np.log10(error_T))+1   #evaluate order of Error
    save_path = 'C:/Users/User/Desktop/Demonstartion/' + name + '/'
    
    fig = plt.figure(I*(Rep+2) + 2)             #plot dots and curve
    fig.set_size_inches(50, 6)
    plt.plot(x, y, '.b')
    plt.xlabel('BJD')
    plt.ylabel('V mmag')
    plt.title('Light curve')
    xx = np.linspace(min(x), max(x), len(x))
    plt.plot(xx, sin(xx, ans_ideal, n), '-r')
    plt.savefig(save_path + name + " light curve.png", dpi = 300)
    plt.show()
    
    return T_ideal, ans_ideal, np.round(T_ideal, order_Error), np.round(error_T, order_Error), ans_ideal

def becoming_perfect_second(I, answ, N, x, y, y_err, name, Numbers, Parametr, n, I_star, Rep):  
    p0 = np.zeros(2*n+1)
    
    for i in range(5):  
        p0[i] = answ[i]
    for i in range(5, 2*n+1):  
        p0[i] = 1
    fun = lambda pp: (y - sin(x, pp, n))/y_err       #core of least squares
    ans = spo.leastsq(fun, p0, full_output=1)
    sigma = np.sqrt(np.sum((y - sin(x, ans[0], n))**2)/len(x))
    error = np.sqrt(np.diag(ans[1]*sigma))
    order_Error = -int(np.log10(error[0]))+1
    
    T_ideal = ans[0][0]  
    Number_of_elements = len(x)
 
    Number_periods = (x - x[0])/T_ideal                   #Phase curve
    X_E = np.zeros(Number_of_elements)
    y_max = y[0]
    for i in range(Number_of_elements):
        if (y[i] > y_max):
            y_max = y[i]
            I_max = i
    for i in range(Number_of_elements):                         ##NEED TO ADD HERE ЗСУВ ПО Х
        X_E[i] = (x[i] - x[0]) - int(Number_periods[i])*T_ideal
    delta = X_E[I_max]
    for i in range(Number_of_elements):                         ##NEED TO ADD HERE ЗСУВ ПО Х
        X_E[i] -= delta
        if (X_E[i] < 0):
            X_E[i] += T_ideal   
            
    save_path = 'C:/Users/User/Desktop/Demonstartion/' + name + '/'
            
    hfont = {'fontname':'Helvetica'}
    fig = plt.figure(3 + I + I_star*(Rep+2))  
    plt.gca().invert_yaxis()            
    fig.set_size_inches(10, 6)
    strin = 'Phase (' + str(np.round(2457000 + x[I_max], 5)) + ' +' + str(np.round(T_ideal, 5)) + ')'
    plt.xlabel(strin, fontsize =14, **hfont)
    plt.ylabel('V, mmag', fontsize = 14, **hfont)
    plt.plot(X_E/T_ideal, y, '.g')
    plt.text(0, (np.min(y)+(1/30)*(np.max(y) - np.min(y))), name, fontsize = 14, **hfont)
    plt.savefig(save_path + name + "phase curve " + str(I) + ".png")
    plt.show()
    
    NName = name + "phase curve" + str(I) + ".txt"
    completeName = os.path.join(save_path, NName) 
    with open(completeName, 'w+') as f:
        for i in range(Number_of_elements):
            f.write(str(X_E[i]) + str(y[i]))

    X_0 = list(x)          
    y_0 = list(y)
    y_0_error = list(y_err)
    k=0
    for i in range(Numbers):
        if ( np.abs( y[i] - sin1(x, ans[0], n, i ) ) > (Parametr*sigma)):
            del X_0[i-k]
            del y_0[i-k]
            del y_0_error[i-k]           
            k+=1       
    x = np.array(X_0)
    y = np.array(y_0)
    y_err = np.array(y_0_error)
    Number_of_elements = len(x)
    return np.round(T_ideal, order_Error), np.round(error[0], order_Error), Number_of_elements, x, y, y_err

"""==========================================="""
"""COMPUTING APPROXIMATE VALUE OF PERIOD"""
"""==========================================="""

def Approximation_T(XxX, YyY, YyY_err, A, I, Rep, TTT_max):
    print('0')
    n_n = 2        #number of additions in Fourie series 
    T_minnn = 0.01
    T_maxxx = TTT_max  #range of periods that can be (not take T_min=0)
    N_N = int(T_maxxx/0.01)       #number of dots in this area
    X_minn = 0      #just for fun(do not change)
    print('1')
    
    def sin(t, T_Tt, pp):           #approximation of function that take x data, period and parametrs and give the approximation function
        x = np.zeros(len(t))     #make array x lenth of x-data and full zero  
        for i in range(n_n):       #additions in Fourie series
            x += pp[2*i]*np.sin(2*np.pi*t/T_Tt*(i+1)+pp[2*i+1])
        return x                 #return tha value of approximation function
    
    def sigma(XxXx, YyYy, YyYy_err, T_Tt):                                         #function to find the sum of squares for each T
        fun = lambda pp: (YyYy - sin(XxXx, T_Tt, pp))/YyYy_err              #core of least squares
        ans = spo.leastsq(fun, p0, full_output=1)
        Sigma = np.sum((YyYy-sin(XxXx, T_Tt, ans[0]))**2)/len(XxX)         #ans[0] - parametrs: amplitudes and phases
        return Sigma 
  
    A = (np.max(YyY) - np.min(YyY))*0.5
    p0 = []                     #start conditions
    p0.append(1)
    p0.append(A)
    for i in range(2*n_n-1):
        p0.append(1)
    YyY_sigma = []                #arrays for finding minimum
    XxX_sigma = []
     
    fig = plt.figure(1+I*(Rep+2))
    fig.set_size_inches(10, 6)
    T_T = np.linspace(T_minnn, T_maxxx, N_N)
    
    for i in T_T:
        XxX_sigma.append(i)
        YyY_sigma.append(sigma(XxX, YyY, YyY_err, i))
        plt.plot(i, sigma(XxX, YyY, YyY_err, i), '.r') 
    
    for i in range(N_N):
        YyY_sigma[i] = YyY_sigma[i]/np.max(YyY_sigma)
    value = False
    
    for i in range(N_N):
        if (YyY_sigma[i] < 0.3) and (not value):
            Index_1 = i
            value = True
        if value:
            if (YyY_sigma[i] > 0.3):
                Index_2 = i
                break

    Y_minn = YyY_sigma[Index_1]
    for i in range(Index_1, Index_2):
        if (YyY_sigma[i] < Y_minn):
            Y_minn = YyY_sigma[i]
            X_minn = XxX_sigma[i]
      
    local_delta = (T_maxxx-T_minnn)/N_N
    order_ld = -int(np.round(np.log10(local_delta)))+1
    #print(np.round(X_minn, order_ld), "+-", np. round(local_delta, order_ld))
    return np.round(X_minn, order_ld), np. round(local_delta, order_ld)

"""==========================================="""
"""CREATING WINDOW AND GENERAL WIDJETS"""
"""==========================================="""
def Manual_work():
    
    """==========================================="""
    """MAIN FUNCTION FOR MANUAL MODE"""
    """==========================================="""
    def do_for_single_star():
        
        start_time = time.time()         #time of begining
    #    name, ftype, Repeats, Part_diagramm, T, dT
        enttime[1].delete(0, len(enttime[1].get()))
        enttime[0].delete(0, len(enttime[0].get()))
        entT.delete(0, len(entT.get()))
        entdT.delete(0, len(entdT.get()))
    
        name = ent_StarName.get()
        ftype = ent_TypeFile.get()
        T = float(init_param[0].get())
        dT = float(init_param[1].get())
        Repeats = int(init_param[2].get())
        Parametr = float(init_param[3].get())
        n = int(init_param[4].get())
        
        try:
        # Create target Directory
            os.mkdir(name)
            print("Directory " , name ,  " Created ") 
        except FileExistsError:
            print("Directory " , name ,  " already exists")        
        
        
        # getting data from file
        x, y, y_err, Number_of_elements, Number_of_elements0 = read_data(name, ftype)
        # setting bounds
        A0 = (max(y)-min(y)) / 2
        bounds = [T+dT, T-dT, 0.8*A0, 1.2*A0, 0, 2*np.pi]
        # calculationg ideal values of the futting parameters  
        T_ideal, ans_ideal, T, ΔT, ans = becoming_perfect(bounds, N, x, y, y_err, name, 10, 100)
        for indicator in range(Repeats+1):
            T, ΔT, Number_of_elements, x, y, y_err = becoming_perfect_second(indicator, ans_ideal, N, x, y, y_err, name, Number_of_elements, Parametr, n, 10, 100)
        t_0 = time.time() - start_time
        enttime[1].insert(0, str(round(t_0)-60*int(t_0/60)))
        enttime[0].insert(0, str(int(t_0/60)))
        entT.insert(0, str(T))
        entdT.insert(0, str(ΔT))
        
    """==========================================="""
    """CLEARING WINDOW"""
    """==========================================="""
    
    def clear_win():
        for i in range(5):
            init_param[i].delete(0, len(init_param[i].get()))
        ent_TypeFile.delete(0, len(ent_TypeFile.get()))
        ent_StarName.delete(0, len(ent_StarName.get()))
        entT.delete(0, len(entT.get()))
        entdT.delete(0, len(entdT.get()))
        for i in range(5):
            enttime[i].delete(0, len(enttime[i].get()))    
        
    window = tnk.Tk()
    bcg_cl = '#9999FF'
    window.title("Period D&P V3.0.1 created 27/05/2019")
    w = 900
    h = 450
    window.geometry(str(w) + 'x' + str(h))
    window.config(bg=bcg_cl)
    window.resizable(width=False, height=False)
    
    lb_head = tnk.Label(window, font = ('Algerian', 19), text = 'Name of star:', bg=bcg_cl)
    
    lb_TypeFile = tnk.Label(window, font = ('Bookman Old Style', 14), text = 'Type of file:', bg=bcg_cl)
    lb_par = tnk.Label(window, font = ('Algerian', 19), text = 'Parameters:', bg=bcg_cl)
    ent_StarName = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    ent_TypeFile = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    
    lb_head.place(x = 20, y = 30)
    lb_TypeFile.place(x = 40, y = 70)
    ent_StarName.place(x = 250, y = 35)
    ent_TypeFile.place(x = 250, y = 72)
    lb_par.place(x = 20, y = 135)
    
    text_init = ['T approximately', 'Error of Τ', 'Repeats', 'Parametr of sigma', 'N additions in Fourie']
    init_param = [tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12) for i in range(5)]
    init_param_lables = [tnk.Label(window, font = ('Century', 14), text = text_init[i], bg=bcg_cl) for i in range(5)]
    for i in range(5):
        init_param_lables[i].place(x = 40, y = 170 + i * 35)
        init_param[i].place(x = 250, y = 175 + i * 35)
        
    lb_Results = tnk.Label(window, font = ('Algerian', 19), text = 'Results:', bg=bcg_cl)
    lbT = tnk.Label(window, font = ('Century', 14), text = 'Period', bg=bcg_cl)
    lbpm = tnk.Label(window, font = ('Century', 14), text = '+-', bg=bcg_cl)
    entT = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 8)
    entdT = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 8)
    time_text = ['Time of calculations', 'min', 's']
    lbtime = [tnk.Label(window, font = ('Century', 14), text = time_text[i], bg=bcg_cl) for i in range(3)]
    enttime = [tnk.Entry(window, font = ('Bookman Old Style', 14), width = 7) for i in range(2)]
    
    lb_Results.place(x= 500, y = 100)
    lbT.place(x = 520, y = 140)
    lbpm.place(x = 700, y = 170)
    entT.place(x = 600, y = 170)
    entdT.place(x = 725, y = 170)
    lbtime[0].place(x = 520, y = 230)
    lbtime[1].place(x = 690, y = 260)
    lbtime[2].place(x = 825, y = 260)
    enttime[0].place(x = 600, y = 260)
    enttime[1].place(x = 735, y = 260)
    
    btn1 = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Calculate', bg = "blue", fg = "white", height = 1, width = 13, command = do_for_single_star)
    btn2 = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Clear', bg = "red", fg = "white", height = 1, width = 13, command = clear_win)  
    btn1.place(x = 200, y = 400)
    btn2.place(x = 20, y = 400)
    window.mainloop()
    
def Automatic_work():
        
    window = tnk.Tk()
    bcg_cl = '#9999FF'
    window.title("Period D&P V3.0.1 created 27/05/2019")
    w = 550
    h = 260
    window.geometry(str(w) + 'x' + str(h))
    window.config(bg=bcg_cl)
    window.resizable(width=False, height=False)
    
    lb_head = tnk.Label(window, font = ('Bookman Old Style', 18), text = 'Task file:', bg=bcg_cl)
    lb_head.place(x = 20, y = 30)
    lb_head = tnk.Label(window, font = ('Bookman Old Style', 15), text = 'Upper evaluation of T:', bg=bcg_cl)
    lb_head.place(x = 20, y = 100)
    
    ent_TaskFile = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    ent_TaskFile.place(x = 70, y = 70)
    ent_Tmax = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 12)
    ent_Tmax.place(x = 70, y = 135)

    
    
    """==========================================="""
    """MAIN FUNCTION FOR AUTOMATIC MODE"""
    """==========================================="""
    
    def automatic_regime():
        enttime[0].delete(0, len(enttime[0].get()))
        enttime[1].delete(0, len(enttime[1].get()))
        ent_progress_1.delete(0, len(ent_progress_1.get()))
        ent_progress_2.delete(0, len(ent_progress_2.get()))
        task_file = ent_TaskFile.get()
        TT_max = float(ent_Tmax.get())
        with open(task_file, 'r') as f:
            s = f.readlines()
        res = ''
        start_time_0 = time.time()
        ent_progress_1.insert(0, '0')
        ent_progress_2.insert(0, str(len(s)))
        for i in range(len(s)):
            try:
                start_time = time.time()
                a = s[i].split()
                name = a[0]
                ftype = a[1]
                Repeats = int(a[2])
                Parametr = float(a[3])
                n = int(a[4])
                
                try:
                # Create target Directory
                    os.mkdir(name)
                    print("Directory " , name ,  " Created ") 
                except FileExistsError:
                    print("Directory " , name ,  " already exists") 
                
                x, y, y_err, Number_of_elements, Number_of_elements0 = read_data(name, ftype)
                A0 = (max(y)-min(y)) / 2
                Tappr, Terr = Approximation_T(x, y, y_err, A0, i, len(s), TT_max)
                bounds = [Tappr - Terr, Tappr + Terr, 0.8*A0, 1.2*A0, 0, 2*np.pi]
                T_ideal, ans_ideal, T, ΔT, ans = becoming_perfect(bounds, N, x, y, y_err, name, i, len(s))
                for indicator in range(Repeats+1):
                    T, ΔT, Number_of_elements, x, y, y_err = becoming_perfect_second(indicator, ans_ideal, N, x, y, y_err, name, Number_of_elements, Parametr, n, i, len(s))               
                t0 = time.time() - start_time
                arr = [name, T, ΔT, np.round(t0, 1)]
                for j in range(4):
                    res += str(arr[j]) + ' '
                res += '\n'
                start_time = t0
                k = int(ent_progress_1.get())
                ent_progress_1.delete(0, len(ent_progress_1.get()))
                ent_progress_1.insert(0, str(k+1))
            except: 
                print("Problem with " + str(i+1) + " star. Please check in manual mode")
                res += 'Problem. Check ' + name + 'manually'
                res += '\n'
        f = open('results.dat', 'w')
        f.writelines(res)
        t_0 = time.time() - start_time_0  
        enttime[1].insert(0, str(round(t_0)-60*int(t_0/60)))
        enttime[0].insert(0, str(int(t_0/60)))
        f.close()
        
    time_text = ['Time of calculations', 'min', 's']
    progress = ['Progress', 'from']
    lbprogress = [tnk.Label(window, font = ('Century', 14), text = progress[i], bg=bcg_cl) for i in range(2)]
    lbtime = [tnk.Label(window, font = ('Century', 14), text = time_text[i], bg=bcg_cl) for i in range(3)]
    enttime = [tnk.Entry(window, font = ('Bookman Old Style', 14), width = 4) for i in range(2)] 
    ent_progress_1 = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 3)
    ent_progress_2 = tnk.Entry(window, font = ('Bookman Old Style', 14), width = 3)
    lbprogress[0].place(x = 290, y = 50)
    lbprogress[1].place(x = 372, y = 80)
    ent_progress_1.place(x = 330, y = 80)
    ent_progress_2.place(x = 420, y = 80)
    lbtime[0].place(x = 290, y = 130)
    lbtime[1].place(x = 377, y = 160)
    lbtime[2].place(x = 480, y = 160)
    enttime[0].place(x = 320, y = 165)
    enttime[1].place(x = 425, y = 165) 
    
    btn = tnk.Button(window, font = ('Bookman Old Style', 14), text = 'Calculate', bg = "blue", fg = "white", height = 1, width = 13, command = automatic_regime)
    btn.place(x = 50, y = 185) 
    window.mainloop()

   
window_0 = tnk.Tk()
bcg_cl = '#9999FF'
window_0.title("Period D&P V3.0.0 created 26/05/2019")
w = 400
h = 100
window_0.geometry(str(w) + 'x' + str(h))
window_0.config(bg=bcg_cl)
window_0.resizable(width=False, height=False)

btn1 = tnk.Button(window_0, font = ('Bookman Old Style', 14), text = 'Manual Work', bg = "blue", fg = "white", height = 1, width = 13, command = Manual_work)
btn2 = tnk.Button(window_0, font = ('Bookman Old Style', 14), text = 'Automatic Work', bg = "green", fg = "white", height = 1, width = 13, command = Automatic_work)
btn1.place(x = 200, y = 30)
btn2.place(x = 20, y = 30)
window_0.mainloop()