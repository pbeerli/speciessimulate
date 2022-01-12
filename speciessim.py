#!/usr/bin/env python
# simulates a tree of samples using the structured coalescence with
# speciation and migration
#
# Haleh Ashki 2014, Tallahassee
# edited by P Beerli 2014
#
######run: python Sim_Spc_Dis_event.py < input_Dist_3n.txt
"""
mainUsage: python sim_species.py #loci input_species.txt seed > treefile
"""

import sys
import random
import math
import time
import dendropy
from dendropy import Tree

numpop=0

# specifies preamblee for the simulation program migdata
# setting a fixed number of sites (1000) per locus
#
def preamble(loci, mutation, numpop, numind, shortpreamble):
    if shortpreamble==False:
        print ("#SN")
        print ("#-1")
        print ("#%li" % numpop)
        print ("#%li" % numind,end=' ')
        for i in range(numpop-1):
            print (numind,end=' ')
        print()
        print ("#%li" % loci,"1000 2.0")
        for i in range(loci):
            print ("# rate among sites for locus 0 (1.000000)")
            print ("#=")
    print ("#$Locus 0")
    print ("#$%.20f" % mutation)

######
# Class Node is used to print the Newick format for the tree.
class Node:
    def __init__(self):
        self.name=-1
        self.left=-1
        self.right=-1
        self.ancestor=-1
        self.branchlength=-1
        self.event=[]


    def debugprint(self):
        print ("name: " + str(self.name))
        print ("left: " + str(self.left))
        print ("right: " + str(self.right))
        print ("blen: " + str(self.branchlength))
        print ("event: " + str(self.event))


    # this prints
    #   (1)    name:branchlength [event][event]...[event]
    #   (2)    :branchlength [event][event]...[event]
    def myprint(self):
        if (self.name != "" ):
            print (self.name,end=' ')
        if (self.branchlength != -1 ):
            print (":%s" %str(self.branchlength),end=' ')
            for i in range(len(self.event)):
                print ("[",end=' ')
                print ("%s" %str(self.event[i]),"]",end=' ')


    def tip(self, name):
        self.name= name


    def interior(self, left, right):
        #print "0mad  to interior"
        #print left, right
        self.name= ""
        self.left= left
        self.right= right
        self.ancestor= -1


class C_Tree(Node):

    def __init__(self, root):
        self.root =root
    def myprinttree(self, p):
        if (p.left != -1):
            print ("(",end='')
            self.myprinttree(p.left)
            print (",",end='')
        if (p.right != -1):
            self.myprinttree(p.right)
            print (")",end='')
        p.myprint()




# a list to store the name of the node and the object_name corresponding to that
Nobject=[]
def find(node_name):
    for j in range(0,len(node)):
        ll=Nobject[j]
        if node_name == ll[0]:
            return ll[1]



# calculates the time u of coalescence and migration
def cal_u(poplist,Mig_pop):
    event_list=[]
    event_prob=[];
    r=random.uniform(0.0,1.0)
    l=math.log(r)
    coal=[]
    for i in range(len(poplist)):
        P=poplist[i];
        ##print "i, Ni", i,P,  N[P-1]
        tempval=float(K[P-1]*(K[P-1]-1))/float(4*N[P-1]);
        coal.append(tempval);
        #print "coal prob", P, i, K[P-1] ,tempval;
        if tempval != 0:
            event_list.append(P);
            event_prob.append(tempval);

    migrate=[];
    for i in range(len(Mig_pop)):
        for j in range(len(Mig_pop)):
            if i != j :
                ##print i, j, Mig_pop[i], Mig_pop[j], float(K[i-1]*float(tableMig[i][j]))
                #tempval=float(K[i-1]*float(tableMig[i][j]));
                tempval=float(K[i]*float(tableMig[i][j]));
                migrate.append(tempval)
                if tempval != 0:
                    tempstr=''.join([str(Mig_pop[i]),':',str(Mig_pop[j])])
                    event_list.append(tempstr);
                    event_prob.append(tempval);


    d=sum(coal)+sum(migrate);

    if d==0 :
        u=10**20;
    else:
        u = -(l/d);

    return [u, event_list,event_prob]


def cal_event(eventlist,eventprob):
    #print "cal_event",eventlist,  eventprob , "\n";
    d=sum(eventprob);
    S=eventprob;
    ss=len(S);
    L=[0]*int(ss);
    for i in range(ss):
        L[i]=S[i]/d;

    for i in range(1,ss):
        L[i]=L[i-1]+L[i];


    r=random.uniform(0.0,1.0)
    for i in range(0,len(L)):
        if L[i]>r:
            return str(eventlist[i])




def coalescent(PNO):

    global interiorNO,K

    PNO=int(PNO)
    tips=[]
    length=len(node);
    #print length;

    for i in range(0,length):
        nn=node[i];
        if nn[0] != '' and  nn[1]==int(PNO)  and nn[3]=='N':   # not coalescent yet
            #print nn[0]
            tips.append(nn[0])

    #print "tips:  ", tips ,"\n"
    # choose two node to coalecsent randomly from the tips list
    #random.seed(10)   ### it can get the fix number or time.time()  need ti import time then
    #print tips, "@",
    ff=random.sample(tips,2);
    #print ff
    ## adding the new coalesnet interior node in list
    tt=int(sample_number)+interiorNO
    interior="i"+str(tt)
    LRnode=(ff[0],ff[1])
    templist=[];
    node.append([interior,int(PNO),0.0,'N',LRnode,templist])

    interiorNO=interiorNO+1;

    ##  set the coalescented node to Y
    length=len(node);
    for i in range(0,length):
        nn=node[i]
        if  nn[0]==ff[0] or nn[0]==ff[1]:
            node[i]=[nn[0],nn[1],nn[2]*mutation,'Y',nn[4],nn[5],(U-nn[2])*mutation]


    ## 2 nodes have coalsecent so th enumber of tips is reduced  population n is the n-1 index in list
    K[PNO-1]=K[PNO-1]-1;




def Migration(u,fromNo,toNo):
    global K
    ##global node
    tips=[]
    length=len(node);
    fromNo=int(fromNo);
    toNo=int(toNo);

    for i in range(0,length):
        nn=node[i]
        if nn[0] != '' and  nn[1]==int(fromNo) and nn[3]=='N':
            tips.append(nn)


    tempval=[];
    templist='';
    length=len(node);
    ff=random.sample(tips,1)
    ##print "------", ff[0][0], fromNo, toNo
    for i in range(0,length):
        ntemp=node[i]
        if  ntemp[0] == ff[0][0]:
                ##print "migration node", ntemp[0],ff[0][2], u
                #print ff[0][2],"@", mutation
            xx =str(ff[0][2]*mutation)
            templist=[''.join(['&M ',str(fromNo-1),' ',str(toNo-1),':',xx])];
            ##templist=[''.join(['&M ',str(fromNo-1),' ',str(toNo-1),':',str(u)])];
            tempval=ntemp[5]+templist;
#                               node[i]=[ntemp[0],toNo,ntemp[2],ntemp[3],ntemp[4],tempval, "@@"]
            node[i]=[ntemp[0],fromNo,ntemp[2],ntemp[3],ntemp[4],tempval, "@@"]
            break;

#       K[fromNo-1]=K[fromNo-1]-1;
#       K[toNo-1]=K[toNo-1]+1;
    K[fromNo-1]=K[fromNo-1]+1;
    K[toNo-1]=K[toNo-1]-1;




def Divergence_point(S,n,fn, SplitNo):   ## S newick format of tree as string  and N= number of populatins
    #S='((A,B),(C,D))';
    t1=Tree.get_from_string(S, schema='newick')


    L=[]
    for nd in t1.postorder_node_iter():
        if nd.taxon is not None:
            L.append(("T",str(nd.taxon)));
        else:
            L.append(("I",nd.label));

    #print L;
    temp_pop=[]
    T=[]
    c=0;
    for i in range(len(L)):
        c=c+1;
        if L[i][0] == "I" and c==3:
            temp_pop.append((L[i-2][1],L[i-1][1]))
            c=0;
        elif  L[i][0] == "I" and c==1:
            T.append(temp_pop);
            temp_pop=[];
            c=0;

    #print T;

    temp_n=[];
    final=[]
    for i in range(len(T)):
        if len(T[i]) > 1:
            for j in range(len(T[i])):
                a=T[i][j][0];
                b=T[i][j][1];
                string=''.join([str(a),':',str(b)])
                if string in divergance:
                    temp_n.append(divergance[string]);
                else:
                    temp_n.append(n+1+len(divergance.keys()));
                    divergance[string]=n+1+len(divergance.keys());


            string=''.join([str(temp_n[0]),':',str(temp_n[1])])
            if string not in divergance:
                divergance[string]=n+1+len(divergance.keys());


            final.append(divergance[string]);

        else: ### means the lenght is 1
            a=T[i][0][0];
            b=T[i][0][1];
            string=''.join([str(a),':',str(b)])
            if string not in divergance:
                divergance[string]=n+1+len(divergance.keys());

            final.append(divergance[string]);



    if len(final) >1 :
        string=''.join([str(final[0]),':',str(final[1])])
        if string not in divergance:
            divergance[string]=n+1+len(divergance.keys());


    #print "divergence point result: ", string, divergance
    temp_split=string.split(":");
    #print "temp_split", temp_split
    if temp_split[0] not in distribution:
        distribution[temp_split[0]]=[divergance[string],fn];
    if temp_split[1] not in distribution:
        distribution[temp_split[1]]=[divergance[string],fn];

    #print "distribution result", distribution
    return divergance[string];



def speciation_time():
    #print" speciation time";
    global node
    length=len(node);
    tips=[];
    Stime=[];
    nodeval=0;
    templist=[];
    for i in range(0,length):
        #print i
        nn=node[i];
        if  nn[3]=='N' and nn[1] != SplitNo:   # not coalescent yet     and not in the last split point
            tips.append(nn[0]);
            pp=nn[1];  ## population Number
            ##print pp, distribution[str(pp)]
            tempval1=distribution[str(pp)];
            tempval=tempval1[1].split(",");  ## it gives mu and sigma values
            randval=random.gauss(float(tempval[0]),float(tempval[1]));
            #print "#@$",tempval1, tempval, randval, U,
            #print U, randval
            count=0
            while randval < U and count<1000:  ##### chekc for < or >
                randval=random.gauss(float(tempval[0]),float(tempval[1]));
                count += 1;
            #print randval, "/#@$"
            if count>=1000:
                break
            Stime.append(abs(randval));

    index=[i[0] for i in sorted(enumerate(Stime), key=lambda x:x[1])]
    if len(index) != 0:   ### means the index is null and there is no node having above condition
        min=Stime[index[0]];
        min=min-U;
        nodeval=tips[index[0]];
    else:
        #print " not above condition", node
        min=0;

    return [min, nodeval];

def speciation(stime,nodeno):

    length=len(node);
    templist='';
    tempval=[];

    for i in range(0,length):
            #print i
        nn=node[i]
        if  nn[0]== nodeno:

            pp=nn[1];
            tempval=distribution[str(pp)];
            newpp=tempval[0];
            #templist=[''.join(['&S ',str(newpp) ,' ',str(pp),':',str(stime)])];
            templist=[''.join(['&S ',str(newpp-1) ,' ',str(pp-1),':',str(nn[2] *mutation )])];

            tempval=nn[5]+templist;
            #print "-----", tempval
            node[i]=[nn[0],newpp,nn[2],nn[3],nn[4],tempval,"##"];

            K[pp-1]=K[pp-1]-1;
            K[newpp-1]=K[newpp-1]+1;

            break;


def population_list():
    pop_event=[];
    length=len(node)
    for i in range(0,length):
        nn=node[i];
        if  nn[3]=='N':   # not coalescent yet
            pp=nn[1];
            if K[pp-1] >= 2:   ## enough sample to coalescent
                pop_event.append(pp);
    d=dict((i,pop_event.count(i)) for i in pop_event)
    #print " population list", pp, d
    return list(d.keys())

def mak_node_list(B):
    global N,K,numpop
    ###making k and n list containg population number and sample numbers
    n=int(B[0]);

    #### making N and K list containing population size and sample size
    c=0;
    for i in range(1, n+1):
        N[c]=int(B[i]);
        c=c+1;

    c=0;
    for i in range(n+1,2*n+1):
        K[c]=int(B[i]);
        c=c+1;

    population_size=sum(N);
    sample_number=sum(K)
    numpop = len(K)

    nodearr=('','')
    Elist=[];  ### containing the list of time and th event happens for this node
    node=[['',0,0.0,'N',nodearr,Elist,0]]*int(sample_number)

    temp=[0]*int(n+1);  ### list of all sample numbers plus having 0 as the first element in order to having intervals between sample numbers
    c=0;
    for i in range(n):
        c=c+1
        temp[c]=K[i];


    c=0;
    for i in range(0,len(temp)-1):
        for j in range(temp[i],temp[i]+temp[i+1]):
            nodename=str(i)+"_"+"n"+str(c+1);
            populationNo=i+1;
            rr=node[c]
            node[c]=[nodename,populationNo,rr[2],rr[3],rr[4],rr[5],0]
            c=c+1;


    return node;


if __name__ == "__main__":
    shortpreamble = False
    if len(sys.argv) == 4:
        loci = int(sys.argv[1])
        myfile = sys.argv[2]
        seedval = sys.argv[3]
    else:
        print ("syntax command #loci|-1 filename seed")
        sys.exit(-1)
    if loci == -1:
        loci = 1
        shortpreamble = True
    start = time.time();
    fff = open(myfile,'rU')
    ####Input data
    # seedval=10;
    random.seed(seedval);
    # random.seed(seedval);
    # test e=2.71828183;


    divergance={}   ###dictionary
    distribution={}   ###dictionary
    c=0
    A=[];
    B=[];
    SplitNo=0;
    for i in fff:
        #print i
        A=i.split(':');
        B.append(A[1])
        if A[0] == 'Divergance Point':
            ###making k and n list containg population number and sample numbers
            n=int(B[0]);

            #### making N and K list containing population size and sample size
            N=[0]*int(n);
            K=[0]*int(n);

            #####making node list
            node=mak_node_list(B)

            population_size=sum(N);
            sample_number=sum(K)
            ####print "+++", node
            interiorNO=1;
            U=0; ##Sum of all time intervals

            B=[];
        if A[0]== "Migration Matrix":
            SplitNo=n;
            for j in range(0,len(B)-3,3):
                SplitNo=SplitNo+1;
                P1=B[j].rstrip('\n');
                P2=B[j+1].rstrip('\n');
                fn=B[j+2].rstrip('\n');

                if len(P1.split(",")) > 1:
                    tt1=int(Divergence_point(str(P1),n,fn,SplitNo))
                else:
                    distribution[P1]=[SplitNo,fn];
                    tt1=P1;
                if len(P2.split(",")) > 1:
                    tt2=int(Divergence_point(str(P2),n,fn,SplitNo))
                else:
                    distribution[P2]=[SplitNo,fn];
                    tt2=P2;


            divergance[''.join([str(tt1),':',str(tt2)])]=SplitNo;
            if str(tt1) not in distribution:
                distribution[str(tt1)]=[SplitNo,fn];
            if str(tt2) not in distribution:
                distribution[str(tt2)]=[SplitNo,fn];

            B=[];
        #PB
        if A[0] == 'mutation':
            mutation = float(A[1])
            #print "@@@ mutation rate = ", mutation,"@@@@@"

        if A[0]=="END":
            tableMig= [ [ 0 for j in range(SplitNo-1) ] for k in range(SplitNo-1) ];
            c=0;
            for j in range(SplitNo-1):
                for k in range(SplitNo-1):
                    tableMig[j][k]=B[c].rstrip('\n');
                    c=c+1;
            break;

    migrate_pop=[x for x in range(1, SplitNo)]
    #### Calculating N(population size) for split points
    ###  I can change this to get thi svalue from user as input ***********
    keylist=list(divergance.keys())
    for i in range(len(keylist)):
        arrval=keylist[i].split(":");
        newN=N[int(arrval[0])-1]+N[int(arrval[1])-1];
        N.append(newN);
        K.append(0);


    ##print "*************", N, K
    counter=0;
    U=0;

    while sum(K) > 1:
        counter=counter+1;
        ##      print "--------------------------------", counter

        pop_list=population_list()

        Result_list=cal_u(pop_list,migrate_pop);
        u=Result_list[0];
        Ev_list=Result_list[1];
        Ev_prob=Result_list[2];

        ##print "u after cal_u", u
        ST_list=speciation_time()


        St=ST_list[0];
        SpcNode=ST_list[1];
        ##print "ST", St
        ucopy=u;

        if ((u> St) and (St != 0)):   ## if that is the speciation then u=U-St
            #               u=St-U;
            u=St;





        U=U+u;

# add new interval to all not_colescent nodes

        length=len(node);
        for i in range(0,length):
            nn=node[i]
            if  nn[3]=='N':
                aa=nn[2]+u
                node[i]=[nn[0],nn[1],aa,nn[3],nn[4],nn[5]]




        if (ucopy> St) and (St != 0):  ##3 it was St instead of u
            speciation(u,SpcNode);
            #print "Speciation", SpcNode, u;

        else: ### u< St
            ##print "events", u,St
            Event=cal_event(Ev_list,Ev_prob)
            #print "@@@@@", Event
            temparr=Event.split(":");
            #print temparr
            if len(temparr) == 1:   ###coalescent event
                #       print "coal";
                coalescent(temparr[0])
            else:   ###migration event
                #       print "migrate", u, U;
                Migration(u,temparr[1],temparr[0]) #ATTENTION
                ##Migration(U,temparr[0],temparr[1])


    Nobject=[]
    for j in range(0,len(node)-1):
        mm=node[j]
        #print "#@mm",mm
        if (mm[5] != []):
            for i in mm[5]:
                temp = i.split(':')
                deltatime = float(temp[1])
                temp = temp[0].split( )
                type = temp[0]
                fromm = int(temp[1])
                to = int(temp[2])
                time = float(mm[6]) + deltatime
        #               print '#@mm',type,fromm,to,time
        Fletter= mm[0][2:3]
        if Fletter == 'n':
            templist=str(mm[1])+'_'+str(mm[0]);
            #print templist;
            Anode=mm[0];
            Anode=Node()
            Anode.branchlength=mm[2]
            Anode.tip(mm[0])
            Anode.event=mm[5]
            bb=(mm[0],Anode)
            Nobject.append(bb)
        else:
            aa=mm[4]
            templist=str(mm[1])+'_'+str(mm[0]);
            Bnode=mm[0];
            Bnode=Node()
            Bnode.branchlength=mm[2]
            Bnode.event=mm[5]
            left=find(aa[0])
            rigth=find(aa[1])
            Bnode.interior(left,rigth)
            bb=(mm[0],Bnode)
            Nobject.append(bb)

    #print "#@mm",node[-1]
    #print U
    #print node                     ##print "\n","------------------------", "\n";
    x=len(node)
    #print "all nodes", node
    mm=node[x-1]
    aa=mm[4]
    Bnode=Node()
    left=find(aa[0])
    rigth=find(aa[1])
    Bnode.interior(left,rigth)
    kk=Bnode

    tree=C_Tree(kk)
    numind = sample_number/numpop
    preamble(loci, mutation, numpop, numind, shortpreamble)
    tree.myprinttree(kk)
    kk=()
    print (";")
