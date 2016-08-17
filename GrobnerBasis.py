from numpy import *
from copy import deepcopy

"""Following is the program to check if a given generator
is a Grobener basis or not for given monomial order."""


def array_to_object(arrayin):
    obj=[]
    for i in list(arrayin):
        obj.append(tuple(i))
    return obj


class Poly(object):
    """takes list of tuples
        each list represent a monomial with first entry as coeficient
        5xy-2x^2 = Poly([(5,1,1),(-2,2,0)])"""
    def __init__(self, poly):
        self.vari=deepcopy(len(poly[0])-1)
        zeropoly=[0]
        zeroterm=(0,)*(self.vari+1)
        zeropoly[0]=zeroterm
        self.zero=zeropoly
        dtype=[("0", float)]
        for i in range(self.vari):
            x=("%d"%(i+1), int)
            dtype.append(x)
        if len(poly)==0:
            self.poly=[]
            return None
        #self.poly=poly 
        d={}
        for i in poly:
            if i[1:] not in d.keys():
                if float(i[0]) != 0.0:
                    d[i[1:]]=i
            else:
                if d[i[1:]][0]+i[0]==0.0:
                    del d[i[1:]]
                else:
                    dummylist=list(d[i[1:]])
                    dummylist[0]=dummylist[0]+i[0]
                    d[i[1:]]=tuple(dummylist)
        unsortpoly=d.values()
        unsortarray=array(unsortpoly, dtype=dtype)
        if len(unsortpoly)==0:
            self.poly=self.zero
            return None
        order=[]
        for i in range(self.vari):
            x="%d"%(i+1)
            order.append(x)
        revsort=sort(unsortarray, order=order)
        reqlist=list(revsort)
        reqlist.reverse()
        self.poly=array_to_object(reqlist)

    def __eq__(self, other):
        return self.poly == other.poly

    def __ne__(self, other):
        return self.poly != other.poly

    def LT(self):                     #note o/p is tuple
        return self.poly[0]

    def leadingterm(self):
        leadingterm=[0]
        leadingterm[0]=self.LT()
        return Poly(leadingterm)

    def multideg(self):
        return self.LT()[1:]

    def isdivisible(self, other):
        for i in range(len(self.multideg())):
            if self.multideg()[i]<other.multideg()[i]:
                return False
        return True

    def __add__(self, other):
        added = self.poly+other.poly
        return Poly(added)

    def __sub__(self, other):
        neglist=[]
        other1=other.poly
        for i in other1:
            dummytuple=(-1*i[0],)+i[1:]
            neglist.append(dummytuple)
        return self.__add__(Poly(neglist))

    def __mul__(self, other):
        if type(other)==int or type(other)==float:
            mullist=[]
            self1=self.poly
            for i in self1:
                dummytuple=(float(other)*i[0],)+i[1:]
                mullist.append(dummytuple)
            return Poly(mullist)
        mullist=[]
        self1=self.poly
        other1=other.poly
        for i in self1:
            for j in other1:
                mulnum=i[0]*j[0]
                muldeg=tuple(array(i[1:])+array(j[1:]))
                multerm=(mulnum,)+muldeg
                mullist.append(multerm)
        return Poly(mullist)

    def monodiv(self, other):
        if type(other)==type(Poly([(0,0,0)])):
            if len(other.poly)==1 and len(self.poly)==1:
                T = self.isdivisible(other)
                if T:
                    a, b = self.poly[0], other.poly[0] 
                    anscoef=float(a[0])/float(b[0])
                    ansdeg=array(a[1:])-array(b[1:])
                    return Poly([(anscoef,)+tuple(ansdeg)])

    def __div__(self, other):
        if type(other)==int or type(other)==float:
            return self.__mul__(1.0/other)             
        s=len(other)   #no. of polynomials
        quotient=[Poly(self.zero)]*s
        remainder=Poly(self.zero)
        dummyself=Poly(self.poly)
        while dummyself.poly != self.zero:
            i=0
            divisionoccurred = False
            while i<s and divisionoccurred == False:
                T = dummyself.isdivisible(other[i])
                if T:
                    quotient[i]=quotient[i]+(dummyself.leadingterm().monodiv(other[i].leadingterm()))
                    dummyself=dummyself-(dummyself.leadingterm().monodiv(other[i].leadingterm()))*other[i]
                    divisionoccurred=True
                else:
                    i=i+1
            if divisionoccurred==False:
                remainder=remainder+dummyself.leadingterm()
                dummyself=dummyself-dummyself.leadingterm()
        return [quotient, remainder]

    def LCM(self, other):
        lead_self=self.LT()
        lead_other=other.LT()
        LCM_term=[1]
        for i in range(len(lead_self)-1):
            LCM_term.append(max(lead_self[i+1], lead_other[i+1]))
        LCM=[tuple(LCM_term)]
        return Poly(LCM)

    def s_poly(self, other):
        LCM = self.LCM(other)
        first = LCM.monodiv(self.leadingterm())
        second = LCM.monodiv(other.leadingterm())
        return first*self-second*other

    def iszero(self):
        return self.poly==self.zero
        
            


"""GB of RNC starts"""
#n is no. of variables-1; varables are x0, x1,...,x_n

def encode_order(term, order=0):  #term is only tuple 
    if order==0:
        return term
    term_list=list(term)
    for i in range(len(order)):
        term_list[i+1]=term[order[i]+1]
    return tuple(term_list)

def minor_helper(n, m, order=0): #gives 2minor object
    if n<=1:
        return "Not Possible"
    if n == 2:
        term_1 = [0]*(m+2)
        term_2 = [0]*(m+2)
        term_1[0], term_1[1], term_1[3] = 1, 1, 1
        term_2[0], term_2[2] = -1, 2
        poly = [encode_order(tuple(term_1), order), encode_order(tuple(term_2), order)]
        return [Poly(poly)]
    minorr=minor_helper(n-1, m, order)
    for i in range(n-2):
        term_1 = [0]*(m+2)
        term_2 = [0]*(m+2)
        term_1[0], term_1[i+1], term_1[n+1] = 1, 1, 1
        term_2[0], term_2[i+2], term_2[n] = -1, 1, 1
        poly = [encode_order(tuple(term_1), order), encode_order(tuple(term_2), order)]
        minorr.append(Poly(poly))
    final_poly_1 = [0]*(m+2)
    final_poly_2 = [0]*(m+2)
    final_poly_1[0], final_poly_1[n-1], final_poly_1[n+1] = 1, 1, 1
    final_poly_2[0], final_poly_2[n] = -1, 2
    final_poly = [encode_order(tuple(final_poly_1), order), encode_order(tuple(final_poly_2), order)]
    minorr.append(Poly(final_poly))
    return minorr

def minor(n, order=0):
    return minor_helper(n, n, order)

def all_s_poly(n, order=0):   #for 2*2 minor problem
    s_poly_list=[]
    minorr=minor(n, order)
    i=0
    while i<n*(n-1)/2-1:
        j=i+1
        while j<n*(n-1)/2:
            s_poly_list.append(minorr[i].s_poly(minorr[j]))
            j=j+1
        i=i+1
    return s_poly_list

def all_s_poly_1(n, p): #p is list of polynomials, n is no of poly
    s_poly_list=[]
    minorr=p
    i=0
    while i<n-1:
        j=i+1
        while j<n:
            s_poly_list.append(minorr[i].s_poly(minorr[j]))
            j=j+1
        i=i+1
    return s_poly_list
    

def isgrobner(n,    order=0):   #for 2*2 minor problem
    return all(map(lambda x: (x/minor(n, order))[1].iszero(), all_s_poly(n, order)))


def isgrobner_1(n, p,    order=0):   #p is list of poly n is no of poly
    return all(map(lambda x: (x/p)[1].iszero(), all_s_poly_1(n, p)))

def f(n, l):
    print "order     :isgrobner"
    for i in l:
        print i,":", isgrobner(n, i)

################Partions##########33
