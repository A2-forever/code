from os.path import exists
from numpy import array
from scipy.special import comb
from scipy.sparse import save_npz
from scipy.sparse import load_npz
from scipy.sparse import lil_matrix as mtr
from scipy.sparse import csr_matrix as mtx

Probability = (1.,.5,.25)
if exists('参数.txt'):
    with open('参数.txt','r') as fp:
        Probability = eval(fp.read())['Probability']


def GetStates():
    StateSpace = []
    surfaces = [(1,1,1),(2,2,2),(1,1,2),(1,1,3),(1,2,2),(2,2,3),(1,2,3)]
    for surface in surfaces:
        for n_1 in range(72-surface.count(1)):
            for n_2 in range(29-surface.count(2)):
                for n_3 in range(2-surface.count(3)):
                    StateSpace.append((surface,(n_1,n_2,n_3)))
    left = [(0,1,3),(0,1,1),(0,0,1),(0,0,3),(0,0,0),(0,2,3),(0,2,2),(0,0,2),(0,1,2)]
    for surface in left:
        StateSpace.append((surface,(0,0,0)))
    StateSpace.sort(key=lambda x:sum(x[1]),reverse=True)
    StateIndex = dict((s,n) for n,s in enumerate(StateSpace))
    return StateSpace,StateIndex

def GetCombs(hidenPool):
    assert sum(hidenPool)>=3
    p_0 = comb(sum(hidenPool),3)
    for n_1 in range(min(3,hidenPool[0])+1):
        p_1 = comb(hidenPool[0],n_1)
        for n_2 in range(min(3,hidenPool[1])+1):
            p_2 = comb(hidenPool[1],n_2)
            for n_3 in range(min(3,hidenPool[2])+1):
                p_3 = comb(hidenPool[2],n_3)
                if n_1 + n_2 + n_3 == 3:
                    p = p_1*p_2*p_3/p_0
                    surface = tuple([*[1]*n_1,*[2]*n_2,*[3]*n_3])
                    pool = tuple(hidenPool-array([n_1,n_2,n_3],dtype=int))
                    state = (surface,pool)
                    yield (state,p)

def Getdata():
    StateSpace,_ = GetStates()
    Surfaces = []
    hidenPool = []
    for state in StateSpace:
        surface = array([state[0].count(i) for i in range(1,4)],dtype=float)
        Surfaces.append(surface)
        hidenPool.append(array(state[1],dtype=float)+surface)
    return array(Surfaces),array(hidenPool)


class ProblemPulse:
    def __init__(self,prefer=2):
        self.prefer = prefer

    def GetSuccessor(self,state):
        arm = 3 if 3 in state[0] else self.prefer
        arm = arm if arm in state[0] else 3 - arm
        surface,pool = list(state[0]),list(state[1])
        if sum(pool) == 0:
            surface.remove(arm)
            surface.append(0)
            nextState = (tuple(sorted(surface)),tuple(pool))
            yield (state,1-Probability[arm-1])
            yield (nextState,Probability[arm-1])
        else:
            surface.remove(arm)
            for i in range(3):
                if pool[i] > 0:
                    p = Probability[arm-1]*pool[i]/sum(pool)
                    nextSurface = sorted([*surface,i+1])
                    nextPool = list(pool)
                    nextPool[i] -= 1
                    nextState = (tuple(nextSurface),tuple(nextPool))
                    yield (nextState,p)
            if Probability[arm-1] < 1:
                for i in range(3):
                    if pool[i] > 0:
                        p = (1-Probability[arm-1])*pool[i]/sum(pool)
                        nextSurface = sorted([*surface,i+1])
                        nextPool = list(pool)
                        nextPool[arm-1] += 1
                        nextPool[i] -= 1
                        nextState = (tuple(nextSurface),tuple(nextPool))
                        yield (nextState,p)
    @staticmethod
    def isGoal(state):
        return 3 not in state[0] and state[1][-1]==0

class ProblemFire(ProblemPulse):
    def GetSuccessor(self,state):
        surface,pool = list(state[0]),state[1]
        if sum(pool) == 0:
            p = 1/(3-surface.count(0))
            for choice in surface[:]:
                if choice != 0:
                    surface.remove(choice)
                    nextSurface = sorted([*surface,0])
                    nextState = (tuple(nextSurface),pool)
                    yield (nextState,p)
                    surface.append(choice)
        else:
            pre = 1/(sum(pool)+3)
            for choice in surface[:]:
                surface.remove(choice)
                for i in range(3):
                    if pool[i] > 0:
                        p = pre*pool[i]/sum(pool)
                        nextSurface = sorted([*surface,i+1])
                        assert len(nextSurface) == 3
                        nextPool = list(pool)
                        nextPool[i] -= 1
                        nextState = (tuple(nextSurface),tuple(nextPool))
                        yield (nextState,p)
                surface.append(choice)
            surface.sort()
            for i in range(3):
                if pool[i] > 0:
                    p = pre*pool[i]
                    nextPool = list(pool)
                    nextPool[i] -= 1
                    nextState = (tuple(surface),tuple(nextPool))
                    yield (nextState,p)

class ProblemFlush(ProblemPulse):
    def GetSuccessor(self,state):
        surface,pool = state
        hidenPool = array([surface.count(i) for i in range(1,4)],dtype=int)+array(pool,dtype=int)
        return GetCombs(hidenPool)
    @staticmethod
    def isGoal(state):
        return (3 not in state[0] and state[1][-1] == 0) or sum(state[1])==0

def GetProbability(problem):
    frontier,index = GetStates()
    Possibles = mtr((len(frontier),len(frontier)))
    for state in frontier:
        if problem.isGoal(state):
            Possibles[index[state],index[state]] = 1.
            continue
        for nextState,p in problem.GetSuccessor(state):
            Possibles[index[state],index[nextState]] += p
    return mtx(Possibles)

def SaveMatrix(new=False):
    if new or not exists('Pulse1.npz'):
        Pulse1 = GetProbability(ProblemPulse(1))
        save_npz('Pulse1.npz',Pulse1)
    if new or not exists('Pulse2.npz'):
        Pulse2 = GetProbability(ProblemPulse(2))
        save_npz('Pulse2.npz', Pulse2)
    if new or not exists('Fire.npz'):
        Fire = GetProbability(ProblemFire())
        save_npz('Fire.npz',Fire)
    if new or not exists('Flush.npz'):
        Flush = GetProbability(ProblemFlush())
        save_npz('Flush.npz',Flush)
    return True

def Estimate(surface,hidenPool,prefer,Pulse,Super,Fire,Flush=0,spac=1):
    assert prefer in [1,2]
    SaveMatrix()
    if sum(surface) < 3:
        if sum(hidenPool) <= 3:
            surface = list(hidenPool)
        else:
            surface = None
    _,indexDict = GetStates()
    hidenPool = array(hidenPool,dtype=int)
    P = mtr((1,len(indexDict)))
    if surface is None:
        for state,p in GetCombs(hidenPool):
            P[0,indexDict[state]] = p
    else:
        assert sum(surface) <= hidenPool.sum()
        pool = tuple(hidenPool-array(surface,dtype=int))
        surface = tuple([*[0]*(3-sum(surface)),*[1]*surface[0],*[2]*surface[1],*[3]*surface[2]])
        state = (surface,pool)
        P[0,indexDict[state]] = 1.
    del indexDict
    P = mtx(P)
    pPulse = load_npz(f'Pulse{prefer}.npz')
    if Flush == 0:
        for _ in range(Pulse):
            P = P.dot(pPulse)
    else:
        spac = max(1,spac)
        pFlush = load_npz('Flush.npz')
        for t in range(0,Pulse,spac):
            for _ in range(min(spac,Pulse-t)):
                P = P.dot(pPulse)
            if Flush > 0 and Pulse-t>=spac:
                P = P.dot(pFlush)
                Flush -= 1
    for _ in range(Super):
        P = P.dot(pPulse)
    pFire = load_npz('Fire.npz')
    for _ in range(Fire):
        P = P.dot(pFire)
    Shown,Pools = Getdata()
    pSurface = P.dot(Shown)
    pGet = hidenPool - P.dot(Pools)
    return pSurface,pGet
