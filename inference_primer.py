import pandas as panda
import scipy.stats as stat
import tkinter as tk
from tkinter import filedialog, MULTIPLE,NORMAL,DISABLED
from tkinter import messagebox
import numpy as np
# import scipy.linalg as la
import itertools as itrt
import pandas as pd

Nsizeforfiting=100.
selectedexp="64" 
automate =False
xlfilename='C:/Users/gboul/pythonProjects/statist/V_2_28S seq all POOLs_edited_no_plummbeus_b_33.xlsx'
selprimer="primer3"
def add_sheet_to_excel(observed, expected, p_value, file_path, sheet_name):
    # Create a DataFrame with the observations, expectations, and p-value
    data = {'Observed': observed, 'Expected': expected, 'P-Value': p_value}
    df = pd.DataFrame(data)

    # Write the DataFrame to the Excel file as a new sheet
    with pd.ExcelWriter(file_path, mode='a', engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name=sheet_name, index=False)
class mystatclass() :
    def __init__(self,sign_level=0.05,onetail=True):
        self.allexp=[]
        self.allobs=[]
        self.allexpP=[]
        self.allobsP=[]
        self.allexpF=[]
        self.allobsF=[]
        self.apprInt=True

        self.con_table = [[]]
        self.sign_level=sign_level
        self.onetail=True
        self.last_test='None'
        self.primer_list=""
        self.simpleaverage=False
        # print ("self.simpleaverage +"+str(self.simpleaverage))

#          Start  implementation methods for reading from excel

    def Scon_table(self,con_table):
        ''' save con_table'''
        self.con_table = con_table

    def conte_table(self):
        ''' exports con_table '''
        return self.con_table 

    def get_primers(self,ic) :
        R = self.R_reads_all[ic]
        primers= list(R.keys())
        return primers

    def form_system(self,ic):
        iexp = self.experiments.index(ic)
        Nsize=0.
        for ip in self.N_mosq_all[ic].keys():
            Nsize+=self.N_mosq_all[ic][ip]
        if Nsize < 1 :
            Nsize=Nsizeforfiting
        R = self.R_reads_all[ic]
        primers= list(R.keys())
        F={}
        for ipr in primers:
            F[ipr] = self.f_fref_est[ipr]

        keys =list(F.keys())
        a= self.formmatrix(keys,R,F)
        print (a)
        A=np.array(a)
#        eigv,eigvec=self.eigenv(A)
        u,s,vh =self.SVD(A) 
        res=self.findzeroeig(s,vh)
        if len(res) > 1 :
            print ("Chech why there are more than one solutions")
        for isol in range(len(res)):
            nres =res[isol][1]/sum(res[isol][1])*Nsize*1.
            print(nres)
            print (np.dot(A,nres))
            inres=[]
            todist=0
            for ire in nres:
                if int(ire)-1 > -1:
                    inres.append(int(ire)-1)
                else :
                    inres.append(int(ire))
                todist+=ire-inres[-1]
            if self.apprInt:
                intsol=self.findintegersol(A,inres,todist)
            else :
                intsol={}

        primedic= {}
        primedicInt= {}

        ike_id = -1
        for ike in keys :
            ike_id +=1 
            primedic[ike] = res[0][1][ike_id]/sum(res[0][1])*Nsize*1.
            if self.apprInt:
                primedicInt[ike] = intsol[ike_id]
        if self.apprInt:
            return primedic,primedicInt
        else :
            return primedic

    def findintegersol(self,A,inres,N):
        # Rgroups=len(inres)
        lx= self.combsolutions(inres,N,Nrange=2)
        bestsol=self.minres(A,lx)
        return lx[bestsol]

    def combsolutions(self,inres,N,Nrange):
        N=int(round(N))
        groups=len(inres)
        allcomb= [ele for ele in itrt.product(range(0, Nrange + 1), repeat = groups)]
        alllist=[]
        for icomn in allcomb :
            if sum(icomn) ==N :
                alllist.append(list(icomn))
        allcomb=[]
        ialllist= 0 
        for alist in alllist:
            allcomb.append(list(np.array(alllist[ialllist])+np.array(inres)))
            ialllist+= 1 
        return allcomb


    def resigual(self,A,x):
        b=np.dot(A,x)
        res=[]
        for ib in b:
            res.append(ib*ib)
        sumres=sum(res)
        return sumres

    def minres(self,A,lx):
        Lres=[]
        for x in lx:
            Lres.append(self.resigual(A,x))
        solumin=Lres.index(min(Lres)) 
        return solumin


    def readxlfileInf(self,filename):
        ''' load read excel file to panda '''
        self.allexp=[]
        self.allobs=[]
        self.allexpP=[]
        self.allobsP=[]
        self.allexpF=[]
        self.allobsF=[]
        self.xlfile=panda.read_excel(filename, sheet_name=None)
        self.experiments = list(self.xlfile.keys())
        self.EHead={}
        self.f_fref_est={}
        for iexp in self.experiments:
            self.EHead[iexp]=self.xlfile[iexp].columns.tolist()
# 1 PRIMNAME	PRIMID	READS	COPIES
        self.Dread_all={}
        self.N_mosq_all={}
        self.Pred={}
        self.R_reads_all={}
        for iexp in self.experiments :
            Dread={}
            for ih in self.EHead[iexp] :
                Dread[ih]=self.xlfile[iexp][ih].tolist()
            try :
                self.Dread_all[iexp]=Dread['PRIMNAME']
            except :
                print (iexp,Dread)
            print ("expe",str(iexp))
        myexperiments=self.experiments
        return myexperiments

    def computefromExcel(self):
        prim_types_all={}
        F_reads_all={}
        self.f_fref_est={}                
        for iexp in self.experiments :
            Dread={}
            for ih in self.EHead[iexp] :
                Dread[ih]=self.xlfile[iexp][ih].tolist()
#            self.Dread_all[iexp]=Dread['PRIMNAME']
            if automate :
                if iexp == selectedexp :
                    Dread["#COPIES"]=Dread["COPIES"]
                    del Dread["COPIES"]
            prim_types={}
            global selectid
            global ref_r  # starting from 1
            print ("add code to select priemer")
            selectid[iexp] = Dread['PRIMNAME'].index(self.primer_list)
            ref_pos[iexp]=selectid[iexp]
            if selectid[iexp] < 0 :
                print ("Reference Prim is not among the ones used in experiment, List starts from 1")

            for ip in Dread['PRIMNAME'] :
                if ip in prim_types.keys() :
                    prim_types[ip]+=1
                else :
                    prim_types[ip]=1
            prim_types_all[iexp]=prim_types

            nd1 = len(Dread['PRIMID'])
            nd2 = len(Dread['READS'])

            if 'COPIES' in Dread.keys() :
                nd3 = len(Dread['COPIES'])
            else :
                nd3 =0
            if '#COPIES' in Dread.keys() :
                nd4 = len(Dread['#COPIES'])
            else :
                nd4 =0

            if nd1 != nd2  :
                print ("Error in data read")
                exit() 
            else :
                if nd3 > 0 :
                    if  nd1 != nd3 :
                        print(nd1)
                        print(nd3)
                        print ("Error in data read")
                        exit() 
            if nd3 > 0 :
                print ("READ the number of Mosquito  ")
                N_mosq =[]
                F_reads=[]
                self.N_mosq_all[iexp]={}
                self.Pred[iexp]=False
                F_reads_all[iexp]={}
                for i in range(nd1): 
                    # can be added to read more Ids starts from 1
                    if Dread['PRIMID'][i] == Dread['PRIMID'][selectid[iexp]] :
                        N_mosq.append([Dread['COPIES'][i],Dread['PRIMNAME'][i]])
                        F_reads.append([Dread['READS'][i],Dread['PRIMNAME'][i]])
                        self.N_mosq_all[iexp][Dread['PRIMNAME'][i]] =Dread['COPIES'][i]
                        F_reads_all[iexp][Dread['PRIMNAME'][i]] =Dread['READS'][i]
            else :
                if nd4 > 0 :
                    print (" Did not READ  number of Mosquito  Can be used to predict ")
                    N_mosq =[]
                    F_reads=[]
                    self.N_mosq_all[iexp]={}
                    self.Pred[iexp]=True
                    F_reads_all[iexp]={}
                    for i in range(nd1): 
                        # can be added to read more Ids starts from 1
                        if Dread['PRIMID'][i] == Dread['PRIMID'][selectid[iexp]] :
                            F_reads.append([Dread['READS'][i],Dread['PRIMNAME'][i]])
                            self.N_mosq_all[iexp][Dread['PRIMNAME'][i]] =Dread['#COPIES'][i]
                            F_reads_all[iexp][Dread['PRIMNAME'][i]] =Dread['READS'][i]
                else :
                    print (" Did not READ  number of Mosquito  Can be used to predict ")
                    N_mosq =[]
                    F_reads=[]
                    self.Pred[iexp]=True
                    self.N_mosq_all[iexp]={}
                    F_reads_all[iexp]={}
                    for i in range(nd1): 
                        # can be added to read more Ids starts from 1
                        if Dread['PRIMID'][i] == Dread['PRIMID'][selectid[iexp]] :
                            F_reads.append([Dread['READS'][i],Dread['PRIMNAME'][i]])
                            F_reads_all[iexp][Dread['PRIMNAME'][i]] =Dread['READS'][i]



            sumREAD =0.
            for ifr in F_reads :
                sumREAD += ifr[0]

            R_reads=[]
            self.R_reads_all[iexp]={}
            for ir in range(len(F_reads)):
                primtype = F_reads[ir][1]
                R_reads.append([F_reads[ir][0]/sumREAD,primtype])
                self.R_reads_all[iexp][primtype]=F_reads[ir][0]/sumREAD
             

            if nd3 > 0 :
                f_fref=[]
                for ir in range(len(F_reads)):
                    f_fref.append(R_reads[ir][0]/float(N_mosq[ir][0]))
                    try :
                        self.f_fref_est[F_reads[ir][1]].append([f_fref[ir],iexp])
                    except :
                        self.f_fref_est[F_reads[ir][1]]=[]
                        self.f_fref_est[F_reads[ir][1]].append([f_fref[ir],iexp])
        myexperiments=self.experiments
        return myexperiments

    def inference_system(self,tf_fref_est,):
        global RExfilename
        try :
            self.R_reads_all
        except :
            self.readxlfileInf(RExfilename)
        
    def solvelin(self,A,B):
        x=np.linalg.solve(A, B)
        return x

    def SVD(self,A):
         u,s,vh=np.linalg.svd(A, full_matrices=True)
         return u,s,vh

    # def eigenv(self,A):
    #     eigvals, eigvecs= la.eig(A)	
    #     return eigvals, eigvecs 

    def avfromdic(self,dic,item):
        templist = []
        for idic in dic[item] :
            templist.append(idic[0])
        temp = np.average(templist)
        temp2 = np.std(templist)
        print (item+ " AV ,STD  :")
        print (temp," " ,temp2)

        return np.average(templist)

    def avRatiofromdic(self,dic,item):
        sumavR = 0.
        isum = 0
        telist=[]
        for idr in  range(len(dic[item])):
            iexp = dic[item][idr][1]
            for jdr in  range(len(dic[self.primer_list])):
                if dic[self.primer_list][jdr][1]==iexp :
                    ref = dic[self.primer_list][jdr][0]
            telist.append(dic[item][idr][0]/ref)
            sumavR +=dic[item][idr][0]/ref
            isum+=1.
        sumavR=sumavR/isum
        return sumavR 

    def formAvf(self,keys,F):
        fk={}
        for ik in keys :
            if self.simpleaverage :
                fk[ik]=self.avfromdic(F,ik)
            else :
                fk[ik]=self.avRatiofromdic(F,ik)
        return fk

    def formmatrix(self,keys,R,F):
        fk=self.formAvf(keys,F)
        A=[]
        i_k =-1
        for ik in keys :
            i_k+=1
            j_k =-1
            A.append([])
            for jk in keys :
                j_k+=1
                if i_k ==j_k :
                    aij = (R[ik]-1.)*fk[jk]
                else :
                    aij =(R[ik])*fk[jk]     #  111111111111
                A[i_k].append(aij)
        return A
                
    def findzeroeig(self,eigv,eigvec):
        res=[]
        for ik in range(len(eigv)):
            if np.allclose(eigv[ik],0.) :
                res.append([eigv[ik],eigvec[ik]])
        return res			
			
    def readxlfileSimple(self,filename):
        ''' load read excel file to panda '''
        self.xlfile=panda.read_excel(filename)
        print (self.xlfile.columns.tolist())
        return self.xlfile.columns.tolist()

    def readcolums(self,slist):
        '''select columes, from  loaded excel file into panda '''
        pdata= panda.DataFrame(self.xlfile, columns= slist)
        rlist =[]
        for i in slist :
            rlist.append(pdata[i].tolist())
        trlist = [[rlist[j][i] for j in range(len(rlist))] for i in range(len(rlist[0]))]            
        return trlist
#          End  implementation methods for reading from excel
#################################################################################################3

#   Symple TK windows
class MyTKwindows():
    def __init__(self,):
        self.mystat = mystatclass()
        self.tkwindows()
        global selectidname
        selectidname = 'Culex_pipiens_PRIMERS'
        global ref_r   # select the primer to use starting from 1
        ref_r = 1

        global selectid
        selectid = {}
        global ref_pos
        ref_pos = {}

        global f_fref_est
        f_fref_est=[]

    def copyclip(self,eve):
        val = eve.widget.get("1.0", 'end-1c') 
        self.root.clipboard_clear() 
        self.root.clipboard_append(val)

    def newwidnow1 (self,iexp,primedic,Ndic,predict=True):  
        if predict or self.chb2.get() : 
            neww1 = tk.Toplevel(self.root)
        table ={}
        iil = -1
        for il in list(Ndic.keys()):
            iil+=1
            leb1 = 'rec'+str(iil)
            table[leb1]={'col1':il,'col2':primedic[il],'col2':Ndic[il],'lebel':leb1}
        # js_table=json.dumps(table)
        if predict or self.chb2.get() : 
            mylebel = tk.Label(neww1,text = "Experement :"+iexp)
            scrollbar = tk.Scrollbar(neww1, orient="vertical")   
#        mylistbox2 = tk.Listbox(neww1,width=100, height=30,yscrollcommand=scrollbar.set )
            mylistbox2 = tk.Text(neww1,width=100, height=30,yscrollcommand=scrollbar.set )
            scrollbar.config(command=mylistbox2.yview)
            scrollbar.pack(side="right", fill="y")
            mylistbox2.bind("<Button-1>", self.copyclip) 
        

            ill= -1
            for il in list(primedic.keys()) :
                ill+=1
                if len(list(Ndic.keys())) >0 :
    #                mylistbox2.insert(ill,str(il)+",  pred :"+ str(primedic[il]) + " ,  used :"+str(Ndic[il])) 
                    mylistbox2.insert(tk.END,str(il)+",  pred :"+ str(primedic[il]) + " ,  used :"+str(Ndic[il])+" \n") 

                else :
    #                mylistbox2.insert(ill,str(il)+",  pred :"+ str(primedic[il]) ) 
                    mylistbox2.insert(tk.END,str(il)+",  pred :"+ str(primedic[il])+" \n" ) 

# error analisis
        obs=[]
        exp=[]
        for ipri in Ndic.keys():
            obs.append(Ndic[ipri]) 
            exp.append(primedic[ipri])
            if predict :
                self.mystat.allexpP.append(Ndic[ipri]) 
                self.mystat.allobsP.append(primedic[ipri])
            else :
                self.mystat.allexpF.append(Ndic[ipri]) 
                self.mystat.allobsF.append(primedic[ipri])

            self.mystat.allexp.append(Ndic[ipri]) 
            self.mystat.allobs.append(primedic[ipri])
        if len(list(Ndic.keys())) >0 and min(exp)>0  :    
            pvalue,chisq,sumerror,N= self.error(obs,exp)
            if automate :
                if iexp == selectedexp :
                   add_sheet_to_excel(obs,exp,pvalue,"excelfileout.xlsx",selectedexp+selprimer)
            if predict or self.chb2.get() : 

                ill+=1
                if predict or self.chb2.get() :             
                    if predict :
                        mylistbox2.insert(tk.END,"Predicted : p value :"+ str(pvalue)+" χ2 value :"+ str(chisq) +" \n") 

                    else :
                        mylistbox2.insert(tk.END,"Used in the Fit  : p value :"+ str(pvalue)+" χ2 value :"+ str(chisq) +" \n") 
        if predict or self.chb2.get() : 
            mylebel.pack(side="top")
            mylistbox2.pack(side="left",fill="both", expand=True)
        


    def newwidnow2 (self):  
        neww2 = tk.Toplevel(self.root)


    def WreadexcelB(self):
        self.WselectList()

    def Wsetavertype(self):
        self.mystat.simpleaverage=self.chb1.get()

    def Wsetshowall(self):    
        self.mystat.showallres=self.chb2.get()

    def WsetapprInt(self):
        self.mystat.apprInt=self.apprInt.get()


    def Winference(self):
        self.WreadListfromExcel()

    def WinferenceEst(self):
        self.mystat.simpleaverage=self.chb1.get()
        self.col_list2=self.mystat.computefromExcel()

        for ic in self.col_list :
            if self.apprInt.get():
                res,resInt = self.mystat.form_system(ic)
            else :
                res = self.mystat.form_system(ic)
            primers = self.mystat.get_primers(ic)
            Nres={}
            for ipri in primers:
                if  ic in self.mystat.N_mosq_all.keys() :
                    if len(self.mystat.N_mosq_all[ic]) >0 :
                        Nres[ipri]=self.mystat.N_mosq_all[ic][ipri]
            if self.apprInt.get() :
                self.newwidnow1(ic,resInt,Nres,predict=self.mystat.Pred[ic])
            else :
                self.newwidnow1(ic,res,Nres,predict=self.mystat.Pred[ic])                    
# obs and exp names are reverced 
        if min(self.mystat.allobs) > 0 :
            pval,chsq,sumerror,N=self.error(self.mystat.allexp,self.mystat.allobs)
            mess= "Pvalue : "+str(pval)+"<(obs-exp)2> : "+str(sumerror/float(N))+"\n"
            print("Pvalue : "+str(pval)+"<(obs-exp)2> : "+str(sumerror/float(N))) 
# obs and exp names are reverced 
        if min(self.mystat.allobsF) > 0 :    
            pval,chsq,sumerror,N1=self.error(self.mystat.allexpF,self.mystat.allobsF)
            if N1 >0 :
                mess+= "Used in Fit  : Pvalue : "+str(pval)+"<(obs-exp)2> : "+str(sumerror/float(N1))+"\n"
                print("Used in Fit  : Pvalue : "+str(pval)+"<(obs-exp)2> : "+str(sumerror/float(N1))) 
# obs and exp names are reverced
        if self.mystat.allobsP!= []:
            if min(self.mystat.allobsP) > 0 :    
                pval,chsq,sumerror,N2=self.error(self.mystat.allexpP,self.mystat.allobsP)
                if N2 >0 :
                    mess+= " Predicted : Pvalue : "+str(pval)+"<(obs-exp)2> : "+str(sumerror/float(N2))+"\n"
                    print(" Predicted : Pvalue : "+str(pval)+"<(obs-exp)2> : "+str(sumerror/float(N2))) 
                    msg=tk.messagebox.showinfo(title=None, message=mess)
        else :
            msg=tk.messagebox.showinfo(title=None, message=mess)

    def error(self,obs,exp,sumerror=0.,N=0):
        if len(obs) != len(exp) :
            print ("len obs != len exp")
            return 0,0,0,0
        chi_squared_stat=0.
        ii = 0
        for i in exp:
            if i==0 :
                del exp[ii]
                del obs[ii]
            ii+=1    
        for i in range(len(obs)):
            chi_squared_stat+=((obs[i]-exp[i])**2)/exp[i]
            sumerror+= (obs[i]-exp[i])**2
        N+=len(obs)
        try : 
            goofnessstat = stat.chisquare(f_obs= obs,f_exp= exp)
            return goofnessstat[1],goofnessstat[0],sumerror,N
        except :
            return False ,False ,False ,False 

    def WreadListfromExcel(self):
        if automate :
            self.xlfilename= xlfilename
        else :
            self.xlfilename= filedialog.askopenfilename()
        self.col_listALL=self.mystat.readxlfileInf(self.xlfilename)
        global RExfilename
        RExfilename=self.xlfilename
        it=-1
        for i in self.col_listALL :
            it=it+1
            stri = i.encode('ascii','replace')
            self.listbox.insert(it,stri)
        self.listbox.bind('<<ListboxSelect>>',lambda event: self.onListclick_event())
        self.listbox.pack()
        self.Tprimer_listALL=[]
        for i in self.col_listALL :
            for ip in self.mystat.Dread_all[i] :
                prname=[]
                for inam in self.Tprimer_listALL : 
                    prname.append(inam[0])
                if ip in prname :
                    posi = prname.index(ip)
                    self.Tprimer_listALL[posi].append(i)
                else :
                    self.Tprimer_listALL.append([ip,i])
        self.primer_listALL=[]
        for ipr in self.Tprimer_listALL :
            if len(ipr) == len(self.col_listALL)+1 :
                self.primer_listALL.append(ipr[0])
        for i in self.primer_listALL :
            it=it+1
            stri = i.encode('ascii','replace')
            self.listboxB.insert(it,stri)
        self.listboxB.bind('<<ListboxSelect>>',lambda event: self.onListclick_eventB())
        self.listboxB.pack(side="right")
        


    def WselectList(self):
        self.xlfilename= filedialog.askopenfilename()
        global RExfilename
        RExfilename=self.xlfilename        
        self.col_listALL=self.mystat.readxlfileSimple(self.xlfilename)
        it=-1
        for i in self.col_listALL :
            it=it+1
            stri = i.encode('ascii','replace')
            self.listbox.insert(it,stri)
        self.listbox.bind('<<ListboxSelect>>',lambda event: self.onListclick_event())
        self.listbox.pack(side=" left")
 
    def onListclick_event(self):
        self.Select_header =[]
        self.col_list=[]
        selection = self.listbox.curselection()
        if len(selection) == 2:
            for isel in selection :
                self.Select_header.append(isel)
            # if self.useall :
            #     self.bButtonexcelE.config(state=NORMAL)
            #     self.bButtonexcel.config(state=NORMAL)
            #     self.bButtonexcel2.config(state=NORMAL)
            for i in self.Select_header:
                self.col_list.append(self.col_listALL[i])
        else :
            for isel in selection :
                self.Select_header.append(isel) 
            for i in self.Select_header:
                self.col_list.append(self.col_listALL[i])                
            
    def onListclick_eventB(self):
        self.Select_headerB =[]
        self.primer_list=[]
        selectionB = self.listboxB.curselection()
        if len(selectionB) == 1:
            for isel in selectionB :
                self.Select_header.append(isel) 
            self.bButtonexcelInf.config(state=DISABLED)
            self.bButtonexcelInf2.config(state=NORMAL)
            self.mystat.primer_list=self.primer_listALL[selectionB[0]]
            if automate : self.WinferenceEst()

            

    def tkwindows(self):
        self.root= tk.Tk()
        self.canv = tk.Canvas(self.root, width = 550, height = 350, bg = 'grey')
        self.lbl = tk.Label(self.root,text = "List of Colums...")
        self.lbl1 = tk.Label(self.root,text = 'Select from Experement')
        self.lbl2 = tk.Label(self.root,text = 'Select Reference Primer')

        scrollbar = tk.Scrollbar(self.root, orient="vertical")   

        self.listbox = tk.Listbox(self.root,selectmode=MULTIPLE,exportselection=False)
        scrollbar.config(command=self.listbox.yview)
        self.listbox.pack(side="left") 
        scrollbar.pack(side="left", fill="y")

        self.checkL1 = tk.Button(self.root,text="select all ",command=self.Wsetselectall)

        scrollbarB = tk.Scrollbar(self.root, orient="vertical")   
        self.listboxB = tk.Listbox(self.root,exportselection=False)
        scrollbarB.config(command=self.listboxB.yview)
        self.listboxB.pack(side="right") 
        scrollbarB.pack(side="right", fill="y")

        self.chb1=tk.IntVar()
        self.checkB = tk.Checkbutton(self.root,text="Don't use reference ratio ",onvalue=1,offvalue=0,variable=self.chb1,command=self.Wsetavertype)
        self.mystat.simpleaverage=self.chb1.get()

        self.chb2=tk.IntVar()
        self.checkB2 = tk.Checkbutton(self.root,text="Show all results ",onvalue=1,offvalue=0,variable=self.chb2,command=self.Wsetshowall)
        self.mystat.showallres=self.chb2.get()


        self.apprInt=tk.IntVar()
        self.checkC = tk.Checkbutton(self.root,text="Approximate Integer solution ",onvalue=1,offvalue=0,variable=self.apprInt,command=self.WsetapprInt)
        self.mystat.apprInt=self.apprInt.get()

        self.lbl1.pack(side="left")
        self.checkL1.pack(side="left")

        self.lbl2.pack(side="right")
        self.checkB.pack()
        self.checkB2.pack()
        self.checkC.pack()
        self.canv.pack()
        self.bButtonexcelInf = tk.Button(text="Read Primer Reads from Excel", command=self.Winference, bg='green', fg='white', font=('helvetica', 12, 'bold'),state=NORMAL)
        self.bButtonexcelInf2 = tk.Button(text="Predict from Primer Reads", command=self.WinferenceEst, bg='green', fg='white', font=('helvetica', 12, 'bold'),state=DISABLED)

        self.canv.create_window(150, 75, window=self.bButtonexcelInf)
        self.canv.create_window(150, 110, window=self.bButtonexcelInf2)

        if automate :
            self.Winference()
            self.Wsetselectall()

    def Wsetselectall(self):
        self.listbox.select_set(0, tk.END)
        self.Select_header =[]
        self.col_list=[]
        selection = self.listbox.curselection()
        if len(selection) == 2:
            for isel in selection :
                self.Select_header.append(isel)
            for i in self.Select_header:
                self.col_list.append(self.col_listALL[i])
        else :
            for isel in selection :
                self.Select_header.append(isel) 
            for i in self.Select_header:
                self.col_list.append(self.col_listALL[i])                

    def looptkwindows(self):    
        self.root.mainloop()


if __name__ == "__main__":
    print ("In the case of fiting the Size of the population is assumed to be ",Nsizeforfiting,"it will affect only the case where integers are requested" )
    mytkwin =MyTKwindows()    
    mytkwin.looptkwindows()
