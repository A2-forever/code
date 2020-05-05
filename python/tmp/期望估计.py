from tkinter import *
from ProbAnalysis import *

def logTraceBack():
    import traceback
    with open('debug.log',"w") as f:
        traceback.print_exc(file=f)
        f.flush()
    print(traceback.format_exc())
    raise

class Surface:
    def __init__(self):
        self.window = Tk()
        self.window.title("铁血捕获期望估计")
        self.window.geometry('395x275')
        self.window.resizable(width=False,height=False)
        self.text()
        self.Pulse = Entry(self.window)
        self.Super = Entry(self.window)
        self.Fire = Entry(self.window)
        self.Flush = Entry(self.window)
        self.Spac = Entry(self.window)

        self.maxPool = [71,28,1]
        self.Pool1 = Spinbox(self.window,from_=0,to=self.maxPool[0])
        self.Pool2 = Spinbox(self.window,from_=0,to=self.maxPool[1])
        self.Pool3 = Spinbox(self.window,from_=0,to=self.maxPool[2])

        self.Show1 = Spinbox(self.window,from_=0,to=3)
        self.Show2 = Spinbox(self.window,from_=0,to=3)
        self.Show3 = Spinbox(self.window,from_=0,to=3)

        self.var1 = StringVar()
        self.var2 = StringVar()
        self.var3 = StringVar()
        self.Exp1 = Label(self.window,textvariable=self.var1,bg='white',relief='sunken')
        self.Exp2 = Label(self.window,textvariable=self.var2,bg='white',relief='sunken')
        self.Exp3 = Label(self.window,textvariable=self.var3,bg='white',relief='sunken')

        self.choice = StringVar()
        self.Choose1 = Radiobutton(self.window,text='一星',variable=self.choice,value='1',command=self.choose)
        self.Choose2 = Radiobutton(self.window,text='二星',variable=self.choice,value='2',command=self.choose)

        self.Compute = Button(self.window,text='计算',command=self.estimate)

        self.Pulse.place(x=75,y=45,width=70,height=30)
        self.Super.place(x=75,y=90,width=70,height=30)
        self.Fire.place(x=75,y=135,width=70,height=30)
        self.Flush.place(x=75,y=180,width=70,height=30)
        self.Spac.place(x=75,y=225,width=70,height=30)

        self.Pool1.place(x=215,y=50,width=40,height=30)
        self.Pool2.place(x=275,y=50,width=40,height=30)
        self.Pool3.place(x=335,y=50,width=40,height=30)

        self.Show1.place(x=215,y=90,width=40,height=30)
        self.Show2.place(x=275,y=90,width=40,height=30)
        self.Show3.place(x=335,y=90,width=40,height=30)

        self.Exp1.place(x=215,y=160,width=100,height=25)
        self.Exp2.place(x=215,y=195,width=100,height=25)
        self.Exp3.place(x=215,y=230,width=100,height=25)

        self.Choose1.place(x=322,y=150)
        self.Choose2.place(x=322,y=170)

        self.Compute.place(x=325,y=200,width=50,height=55)

        self.Pulse.insert(0,'56')
        self.Super.insert(0,'4')
        self.Fire.insert(0,'4')
        self.Flush.insert(0,'0')
        self.Spac.insert(0,'0')
        self.Pool1.delete(0,END)
        self.Pool1.insert(0,'71')
        self.Pool2.delete(0,END)
        self.Pool2.insert(0,'28')
        self.Pool3.delete(0,END)
        self.Pool3.insert(0,'1')
        self.choice.set('2')
        self.prefer = 2

    def text(self):
        Label(self.window,text='资源与参数',font=("",12)).place(x=40,y=10)
        Label(self.window,text='电子脉冲',font=("",10)).place(x=15,y=50)
        Label(self.window,text='超导脉冲',font=("",10)).place(x=15,y=95)
        Label(self.window,text='火神重工',font=("",10)).place(x=15,y=140)
        Label(self.window,text='刷新次数',font=("",10)).place(x=15,y=185)
        Label(self.window,text='刷新间隔',font=("",10)).place(x=15,y=230)
        Label(self.window,text='铁血数目',font=("",12)).place(x=265,y=10)
        Label(self.window, text='一星',font=("",10)).place(x=218,y=30)
        Label(self.window, text='二星',font=("",10)).place(x=278,y=30)
        Label(self.window, text='三星',font=("",10)).place(x=338,y=30)
        Label(self.window,text='总数',font=("",10)).place(x=175,y=55)
        Label(self.window,text='显示',font=("",10)).place(x=175,y=100)
        Label(self.window,text='捕获期望',font=("",12)).place(x=228,y=130)
        Label(self.window,text='一星',font=("",10)).place(x=175,y=165)
        Label(self.window,text='二星',font=("",10)).place(x=175,y=200)
        Label(self.window,text='三星',font=("",10)).place(x=175,y=235)

    def choose(self):
        self.prefer = 2 if self.choice.get() is '2' else 1

    def estimate(self):
        Variates = []
        for variate in [self.Pulse,self.Super,self.Fire,self.Flush,self.Spac]:
            try:
                value = int(variate.get())
            except ValueError:
                value = 0
            variate.delete(0,END)
            variate.insert(0,str(value))
            Variates.append(value)
        hidenPool = []
        Shown = []
        Couples = [(self.Pool1,self.Show1),(self.Pool2,self.Show2),(self.Pool3,self.Show3)]
        for i in range(3):
            pool,show = Couples[i]
            try:
                value1 = max(0,min(self.maxPool[i],int(pool.get())))
            except ValueError:
                value1 = self.maxPool[i]
            pool.delete(0,END)
            pool.insert(0,str(value1))
            hidenPool.append(value1)
            value = min(3-sum(Shown),value1)
            try:
                value2 = max(0,min(value,int(show.get())))
            except ValueError:
                value2 = value
            show.delete(0,END)
            show.insert(0,str(value2))
            Shown.append(value2)
        _,expaction = Estimate(Shown,hidenPool,self.prefer,*Variates)
        self.var1.set(str(round(expaction[0,0],10)))
        self.var2.set(str(round(expaction[0,1],10)))
        self.var3.set(str(round(expaction[0,2],11)))


if __name__ == '__main__':
    try:
        Surface().window.mainloop()
    except:
        logTraceBack()