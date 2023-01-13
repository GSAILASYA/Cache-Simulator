import math

input_list=input().split()
BLOCKSIZE=int(input_list[0])
L1_SIZE=int(input_list[1])
L1_ASSOC=int(input_list[2])
L2_SIZE=int(input_list[3])
L2_ASSOC=int(input_list[4])
REPLACEMENT_POLICY=int(input_list[5])
INCLUSION_PROPERTY=int(input_list[6])
trace_file=input_list[7]

inclusive_bit_value=[0]
#calculations for L1
L1_set_size = (L1_SIZE)//(BLOCKSIZE*L1_ASSOC)
L1_exponent_set_size = int(math.log2(L1_set_size))
L1_tag_limit = 32-L1_exponent_set_size-int(math.log2(BLOCKSIZE))
L1_parameters={'L1_reads':0,'L1_writes':0,'L1_readmiss':0,'L1_writemiss':0,'L1_hits':0,'L1_writebacks':0}

L1={}
for i in range(L1_set_size):
    L1[i]=[]

LRU_L1={}
for i in range(L1_set_size):
    LRU_L1[i]={}

FIFO_L1={}
for i in range(L1_set_size):
    FIFO_L1[i]={}

temp2=[]
mid_list=[]

with open(trace_file) as file_open:
    lines = file_open.readlines()

for i in range(len(lines)):
    temp2.append(lines[i].strip('\n').split(' '))
    if len(temp2[i][1])<8:
        number_of_zeroes=8-len(temp2[i][1])
        temp2[i][1]=str(('0'*number_of_zeroes)+str(temp2[i][1]))
        #binary conversion -> 32
    temp2[i][1]=bin(int(temp2[i][1],16))[2:].zfill(32)
    L1_tag_value = temp2[i][1][0:L1_tag_limit+L1_exponent_set_size]
    mid_list.append(L1_tag_value)



Optimal_L1={}

for i in range(len(mid_list)):
    if mid_list[i] in Optimal_L1:
        Optimal_L1[mid_list[i]].append(i)
    else:
        Optimal_L1[mid_list[i]]=[i]

for tag_keys in Optimal_L1.keys():
    Optimal_L1[tag_keys].append(9999999999)


L1_dirty_dict={}
for i in range(L1_set_size):
    L1_dirty_dict[i]={}

L1_optimal_dict={}
for i in range(L1_set_size):
    L1_optimal_dict[i]={}

temp=[]


def LRU_L1_create(L1_set_value): 
    evict_value =  min(LRU_L1[L1_set_value].values())
    victim_tag = [s for s in LRU_L1[L1_set_value] if LRU_L1[L1_set_value][s]==evict_value]
    return victim_tag[0]

def FIFO_L1_create(L1_set_value):
    evict_value =  min(FIFO_L1[L1_set_value].values())
    victim_tag = [s for s in FIFO_L1[L1_set_value] if FIFO_L1[L1_set_value][s]==evict_value]
    return victim_tag[0]

def Optimal_L1_create(L1_set_value):
    evict_value =  max(L1_optimal_dict[L1_set_value].values())
    victim_tag = [s for s in L1[L1_set_value] if L1_optimal_dict[L1_set_value][s]==evict_value]
    return victim_tag[0]

L2_parameters={'L2_reads':0,'L2_writes':0,'L2_readmiss':0,'L2_writemiss':0,'L2_hits':0,'L2_writebacks':0}
L2_tag_info=[]

#calculations for L2
if (L2_SIZE !=0 and L2_ASSOC!=0):
    L2_set_size = (L2_SIZE)//(BLOCKSIZE*L2_ASSOC)
    L2_exponent_set_size = int(math.log2(L2_set_size))
    L2_tag_limit = 32-L2_exponent_set_size-int(math.log2(BLOCKSIZE))
    L2_tag_info.append(L2_set_size)
    L2_tag_info.append(L2_exponent_set_size)
    L2_tag_info.append(L2_tag_limit)


    L2={}
    for i in range(L2_set_size):
        L2[i]=[]

    LRU_L2={}
    for i in range(L2_set_size):
        LRU_L2[i]={}

    FIFO_L2={}
    for i in range(L2_set_size):
        FIFO_L2[i]={}

    L2_dirty_dict={}
    for i in range(L2_set_size):
        L2_dirty_dict[i]={}

    def LRU_L2_create(L2_set_value): 
        evict_value =  min(LRU_L2[L2_set_value].values())
        victim_tag = [s for s in LRU_L2[L2_set_value] if LRU_L2[L2_set_value][s]==evict_value]
        return victim_tag[0]

    def FIFO_L2_create(L2_set_value):
        evict_value =  min(FIFO_L2[L2_set_value].values())
        victim_tag = [s for s in FIFO_L2[L2_set_value] if FIFO_L2[L2_set_value][s]==evict_value]
        return victim_tag[0]
      
def victim_tag_conversion_L1(L1_set_index, victim_tag):
    intermediate=[]
    L1_set_index=bin(L1_set_index)[2:]
    if len(L1_set_index)<L1_exponent_set_size:
        L1_set_index=str(('0'*(L1_exponent_set_size-len(L1_set_index)))+str(L1_set_index))
    victim_tag=victim_tag + L1_set_index
    L2_tag_value = victim_tag[0:L2_tag_info[2]]
    L2_set_value = int(victim_tag[int(L2_tag_info[2]):int(L2_tag_info[2]+L2_tag_info[1])],2)
    intermediate.append(L2_set_value)
    intermediate.append(L2_tag_value)
    return intermediate


def victim_tag_conversion_L2(L2_set_index, victim_tag):
    intermediate=[]
    L2_set_index=bin(L2_set_index)[2:]
    if len(L2_set_index)<L2_tag_info[1]:
        L2_set_index=str(('0'*(int(L2_tag_info[1])-len(L2_set_index)))+str(L2_set_index))
    victim_tag=victim_tag + L2_set_index
    L1_tag_value = victim_tag[0:L1_tag_limit]
    L1_set_value = int(victim_tag[L1_tag_limit:L1_tag_limit+L1_exponent_set_size],2)
    intermediate.append(L1_set_value)
    intermediate.append(L1_tag_value)
    return intermediate


def L1_calling(L1_set, L1_tag):
    if L1_tag in L1[L1_set]:
        temporary_variable = L1[L1_set].index(L1_tag)
        if L1_dirty_dict[L1_set][L1_tag]=='D':
            inclusive_bit_value[0]+=1
        L1[L1_set].pop(temporary_variable)
        if REPLACEMENT_POLICY == 0:
            del LRU_L1[L1_set][L1_tag]
        if REPLACEMENT_POLICY == 1:
            del FIFO_L1[L1_set][L1_tag]
        del L1_dirty_dict[L1_set][L1_tag]
    


def L2_calling(L2_set_value,L2_tag_value,operation):
    if operation=='r':
        L2_parameters['L2_reads']+=1
        if L2_tag_value in L2[L2_set_value]:
            L2_parameters['L2_hits']+=1
            if (REPLACEMENT_POLICY==0):
                LRU_L2[L2_set_value][L2_tag_value]=max(LRU_L2[L2_set_value].values())+1
    
        else:    
            L2_parameters['L2_readmiss']+=1
            if (len(L2[L2_set_value])<L2_ASSOC):
                L2[L2_set_value].append(L2_tag_value)
                if (REPLACEMENT_POLICY==0):
                    if len(LRU_L2[L2_set_value]) != 0:
                        LRU_L2[L2_set_value][L2_tag_value]=max(LRU_L2[L2_set_value].values())+1
                    else:
                        LRU_L2[L2_set_value][L2_tag_value]=0
                
                elif (REPLACEMENT_POLICY==1):
                    if len(FIFO_L2[L2_set_value]) != 0:
                        FIFO_L2[L2_set_value][L2_tag_value]=max(FIFO_L2[L2_set_value].values())+1
                    else:
                        FIFO_L2[L2_set_value][L2_tag_value]=0

                L2_dirty_dict[L2_set_value][L2_tag_value]='NA'
            else:
                if (REPLACEMENT_POLICY==0):
                    LRU_L2[L2_set_value][L2_tag_value]=max(LRU_L2[L2_set_value].values())+1
                    victim_tag = LRU_L2_create(L2_set_value)
                    if INCLUSION_PROPERTY == 1:
                        L1_victim_address = victim_tag_conversion_L2(L2_set_value, victim_tag)
                        L1_calling(L1_victim_address[0], L1_victim_address[1])
                    for i in range(len(L2[L2_set_value])):
                        if L2[L2_set_value][i] == victim_tag:
                            L2[L2_set_value][i] = L2_tag_value                     
                            break
                
                    del LRU_L2[L2_set_value][victim_tag]
                    if L2_dirty_dict[L2_set_value][victim_tag]=='D':
                        L2_parameters['L2_writebacks']+=1
                    
                    L2_dirty_dict[L2_set_value][L2_tag_value]='NA'
                    del L2_dirty_dict[L2_set_value][victim_tag]
                      
                if (REPLACEMENT_POLICY==1):
                    FIFO_L2[L2_set_value][L2_tag_value]=max(FIFO_L2[L2_set_value].values())+1
                    victim_tag = FIFO_L2_create(L2_set_value)
                    if(INCLUSION_PROPERTY==1):
                        L1_victim_address = victim_tag_conversion_L2(L2_set_value, victim_tag)
                        L1_calling(L1_victim_address[0], L1_victim_address[1])
                    for i in range(len(L2[L2_set_value])):
                        if L2[L2_set_value][i] == victim_tag:
                            L2[L2_set_value][i] = L2_tag_value                      
                            break
                    del FIFO_L2[L2_set_value][victim_tag]
                    if L2_dirty_dict[L2_set_value][victim_tag]=='D':
                        L2_parameters['L2_writebacks']+=1
                    L2_dirty_dict[L2_set_value][L2_tag_value]='NA'
                    del L2_dirty_dict[L2_set_value][victim_tag]

    else:
        L2_parameters['L2_writes']+=1
        if L2_tag_value in L2[L2_set_value]:
            L2_parameters['L2_hits']+=1
            if (REPLACEMENT_POLICY==0):
                LRU_L2[L2_set_value][L2_tag_value]=max(LRU_L2[L2_set_value].values())+1   
            L2_dirty_dict[L2_set_value][L2_tag_value]='D'   
        else:
            L2_parameters['L2_writemiss']+=1
            if (len(L2[L2_set_value])<L2_ASSOC):
                L2[L2_set_value].append(L2_tag_value)
                if (REPLACEMENT_POLICY==0):
                    if len(LRU_L2[L2_set_value]) != 0:
                        LRU_L2[L2_set_value][L2_tag_value]=max(LRU_L2[L2_set_value].values())+1
                    else:
                        LRU_L2[L2_set_value][L2_tag_value]=0
                
                elif (REPLACEMENT_POLICY==1):
                    if len(FIFO_L2[L2_set_value]) != 0:
                        FIFO_L2[L2_set_value][L2_tag_value]=max(FIFO_L2[L2_set_value].values())+1
                    else:
                        FIFO_L2[L2_set_value][L2_tag_value]=0

                L2_dirty_dict[L2_set_value][L2_tag_value]='D'
            
            else:
                if (REPLACEMENT_POLICY==0):
                    LRU_L2[L2_set_value][L2_tag_value]=max(LRU_L2[L2_set_value].values())+1
                    victim_tag = LRU_L2_create(L2_set_value)
                    if(INCLUSION_PROPERTY==1):
                        L1_victim_address = victim_tag_conversion_L2(L2_set_value, victim_tag)
                        L1_calling(L1_victim_address[0], L1_victim_address[1])

                    for i in range(len(L2[L2_set_value])):
                        if L2[L2_set_value][i] == victim_tag:
                            L2[L2_set_value][i] = L2_tag_value
                            break
                    del LRU_L2[L2_set_value][victim_tag]
                    if L2_dirty_dict[L2_set_value][victim_tag]=='D':
                        L2_parameters['L2_writebacks']+=1
                    del L2_dirty_dict[L2_set_value][victim_tag]
                    L2_dirty_dict[L2_set_value][L2_tag_value] = "D"  

                if (REPLACEMENT_POLICY==1):
                    FIFO_L2[L2_set_value][L2_tag_value]=max(FIFO_L2[L2_set_value].values())+1
                    victim_tag = FIFO_L2_create(L2_set_value)
                    if(INCLUSION_PROPERTY==1):
                        L1_victim_address = victim_tag_conversion_L2(L2_set_value, victim_tag)
                        L1_calling(L1_victim_address[0], L1_victim_address[1])
                    for i in range(len(L2[L2_set_value])):
                        if L2[L2_set_value][i] == victim_tag:
                            L2[L2_set_value][i] = L2_tag_value
                            break
                    del FIFO_L2[L2_set_value][victim_tag]
                    if L2_dirty_dict[L2_set_value][victim_tag]=='D':
                        L2_parameters['L2_writebacks']+=1
                    del L2_dirty_dict[L2_set_value][victim_tag]
                    L2_dirty_dict[L2_set_value][L2_tag_value] = "D"



for i in range(len(lines)):
    temp.append(lines[i].strip('\n').split(' '))
    if len(temp[i][1])<8:
        number_of_zeroes=8-len(temp[i][1])
        temp[i][1]=str(('0'*number_of_zeroes)+str(temp[i][1]))
        #binary conversion -> 32
    temp[i][1]=bin(int(temp[i][1],16))[2:].zfill(32)
    L1_tag_value = temp[i][1][0:L1_tag_limit]
    if len(temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size])==0:
        L1_set_value=0
    else:
        L1_set_value = int(temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size],2)
    if (L2_SIZE !=0 and L2_ASSOC!=0):
        L2_tag_value = temp[i][1][0:L2_tag_limit]
        L2_set_value = int(temp[i][1][L2_tag_limit:L2_tag_limit+L2_exponent_set_size],2)

    if temp[i][0]=='r':
        L1_parameters['L1_reads']+=1
        if L1_tag_value in L1[L1_set_value]:
            L1_parameters['L1_hits']+=1
            if (REPLACEMENT_POLICY==0):
                LRU_L1[L1_set_value][L1_tag_value]=max(LRU_L1[L1_set_value].values())+1
            if (REPLACEMENT_POLICY==2):
                Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][1:]
                L1_optimal_dict[L1_set_value][L1_tag_value]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][0]
            
        else:    
            L1_parameters['L1_readmiss']+=1
            if (len(L1[L1_set_value])<L1_ASSOC):
                L1[L1_set_value].append(L1_tag_value)
                if (REPLACEMENT_POLICY==0):
                    if len(LRU_L1[L1_set_value]) != 0:
                        LRU_L1[L1_set_value][L1_tag_value]=max(LRU_L1[L1_set_value].values())+1
                    else:
                        LRU_L1[L1_set_value][L1_tag_value]=0
                
                elif (REPLACEMENT_POLICY==1):
                    if len(FIFO_L1[L1_set_value]) != 0:
                        FIFO_L1[L1_set_value][L1_tag_value]=max(FIFO_L1[L1_set_value].values())+1
                    else:
                        FIFO_L1[L1_set_value][L1_tag_value]=0
                
                elif (REPLACEMENT_POLICY==2):
                    Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][1:]
                    L1_optimal_dict[L1_set_value][L1_tag_value]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][0]
                    

                L1_dirty_dict[L1_set_value][L1_tag_value]='NA'
                #print(L1_dirty_dict[L1_set_value])
            else:
                if (REPLACEMENT_POLICY==0):
                    LRU_L1[L1_set_value][L1_tag_value]=max(LRU_L1[L1_set_value].values())+1
                    victim_tag = LRU_L1_create(L1_set_value)
                    for i in range(len(L1[L1_set_value])):
                        if L1[L1_set_value][i] == victim_tag:
                            L1[L1_set_value][i] = L1_tag_value                  
                            break
                    del LRU_L1[L1_set_value][victim_tag]
                    if L1_dirty_dict[L1_set_value][victim_tag]=='D':
                        L1_parameters['L1_writebacks']+=1
                        if (L2_SIZE!=0 and L2_ASSOC!=0):
                            L2_victim_address = victim_tag_conversion_L1(L1_set_value, victim_tag)
                            L2_calling(L2_victim_address[0], L2_victim_address[1],'w')

                    L1_dirty_dict[L1_set_value][L1_tag_value]='NA'
                    del L1_dirty_dict[L1_set_value][victim_tag]
                      
                if (REPLACEMENT_POLICY==1):
                    FIFO_L1[L1_set_value][L1_tag_value]=max(FIFO_L1[L1_set_value].values())+1
                    victim_tag = FIFO_L1_create(L1_set_value)
                    for i in range(len(L1[L1_set_value])):
                        if L1[L1_set_value][i] == victim_tag:
                            L1[L1_set_value][i] = L1_tag_value                      
                            break
                    del FIFO_L1[L1_set_value][victim_tag]
                    if L1_dirty_dict[L1_set_value][victim_tag]=='D':
                        L1_parameters['L1_writebacks']+=1
                        if (L2_SIZE!=0 and L2_ASSOC!=0):
                            L2_victim_address = victim_tag_conversion_L1(L1_set_value, victim_tag)
                            L2_calling(L2_victim_address[0], L2_victim_address[1],'w')
                    L1_dirty_dict[L1_set_value][L1_tag_value]='NA'
                    del L1_dirty_dict[L1_set_value][victim_tag]
                
                if (REPLACEMENT_POLICY==2):
                    victim_tag = Optimal_L1_create(L1_set_value)
                    for j in range(len(L1[L1_set_value])):
                        if L1[L1_set_value][j] == victim_tag:
                            L1[L1_set_value][j] = L1_tag_value                      
                            break
                    del L1_optimal_dict[L1_set_value][victim_tag]
                    if L1_dirty_dict[L1_set_value][victim_tag]=='D':
                        L1_parameters['L1_writebacks']+=1
                        if (L2_SIZE!=0 and L2_ASSOC!=0):
                            L2_victim_address = victim_tag_conversion_L1(L1_set_value, victim_tag)
                            L2_calling(L2_victim_address[0], L2_victim_address[1],'w')
                    
                    del L1_dirty_dict[L1_set_value][victim_tag]
                    Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][1:]
                    L1_optimal_dict[L1_set_value][L1_tag_value]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][0]
                    L1_dirty_dict[L1_set_value][L1_tag_value]='NA'
            if (L2_SIZE!=0 and L2_ASSOC!=0):
                L2_calling(L2_set_value,L2_tag_value,'r')               


    else:
        L1_parameters['L1_writes']+=1
        if L1_tag_value in L1[L1_set_value]:
            L1_parameters['L1_hits']+=1
            if (REPLACEMENT_POLICY==0):
                LRU_L1[L1_set_value][L1_tag_value]=max(LRU_L1[L1_set_value].values())+1   
            if (REPLACEMENT_POLICY==2):
                Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][1:]
                L1_optimal_dict[L1_set_value][L1_tag_value]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][0]
            L1_dirty_dict[L1_set_value][L1_tag_value]='D'   
        else:
            L1_parameters['L1_writemiss']+=1
            if (len(L1[L1_set_value])<L1_ASSOC):
                L1[L1_set_value].append(L1_tag_value)
                if (REPLACEMENT_POLICY==0):
                    if len(LRU_L1[L1_set_value]) != 0:
                        LRU_L1[L1_set_value][L1_tag_value]=max(LRU_L1[L1_set_value].values())+1
                    else:
                        LRU_L1[L1_set_value][L1_tag_value]=0
                
                elif (REPLACEMENT_POLICY==1):
                    if len(FIFO_L1[L1_set_value]) != 0:
                        FIFO_L1[L1_set_value][L1_tag_value]=max(FIFO_L1[L1_set_value].values())+1
                    else:
                        FIFO_L1[L1_set_value][L1_tag_value]=0
                
                elif (REPLACEMENT_POLICY==2):
                    Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][1:]
                    L1_optimal_dict[L1_set_value][L1_tag_value]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][0]

                L1_dirty_dict[L1_set_value][L1_tag_value]='D'
            
            else:
                if (REPLACEMENT_POLICY==0):
                    LRU_L1[L1_set_value][L1_tag_value]=max(LRU_L1[L1_set_value].values())+1
                    victim_tag = LRU_L1_create(L1_set_value)
                    for i in range(len(L1[L1_set_value])):
                        if L1[L1_set_value][i] == victim_tag:
                            L1[L1_set_value][i] = L1_tag_value
                            break
                    del LRU_L1[L1_set_value][victim_tag]
                    if L1_dirty_dict[L1_set_value][victim_tag]=='D':
                        L1_parameters['L1_writebacks']+=1
                        if (L2_SIZE!=0 and L2_ASSOC!=0):
                            L2_victim_address = victim_tag_conversion_L1(L1_set_value, victim_tag)
                            L2_calling(L2_victim_address[0], L2_victim_address[1],'w')
                    del L1_dirty_dict[L1_set_value][victim_tag]
                    L1_dirty_dict[L1_set_value][L1_tag_value] = "D"  

                if (REPLACEMENT_POLICY==1):
                    FIFO_L1[L1_set_value][L1_tag_value]=max(FIFO_L1[L1_set_value].values())+1
                    victim_tag = FIFO_L1_create(L1_set_value)
                    for i in range(len(L1[L1_set_value])):
                        if L1[L1_set_value][i] == victim_tag:
                            L1[L1_set_value][i] = L1_tag_value
                            break
                    del FIFO_L1[L1_set_value][victim_tag]
                    if L1_dirty_dict[L1_set_value][victim_tag]=='D':
                        L1_parameters['L1_writebacks']+=1
                        if (L2_SIZE!=0 and L2_ASSOC!=0):
                            L2_victim_address = victim_tag_conversion_L1(L1_set_value, victim_tag)
                            L2_calling(L2_victim_address[0], L2_victim_address[1],'w')
                    del L1_dirty_dict[L1_set_value][victim_tag]
                    L1_dirty_dict[L1_set_value][L1_tag_value] = "D"
                
                if (REPLACEMENT_POLICY==2):
                    
                    victim_tag = Optimal_L1_create(L1_set_value)
                    for j in range(len(L1[L1_set_value])):
                        if L1[L1_set_value][j] == victim_tag:
                            L1[L1_set_value][j] = L1_tag_value
                            break
                    del L1_optimal_dict[L1_set_value][victim_tag]
                    if L1_dirty_dict[L1_set_value][victim_tag]=='D':
                        L1_parameters['L1_writebacks']+=1
                        if (L2_SIZE!=0 and L2_ASSOC!=0):
                            L2_victim_address = victim_tag_conversion_L1(L1_set_value, victim_tag)
                            L2_calling(L2_victim_address[0], L2_victim_address[1],'w')
                    del L1_dirty_dict[L1_set_value][victim_tag]
                   
                    Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][1:]
                    L1_optimal_dict[L1_set_value][L1_tag_value]=Optimal_L1[L1_tag_value+temp[i][1][L1_tag_limit:L1_tag_limit+L1_exponent_set_size]][0]
                    L1_dirty_dict[L1_set_value][L1_tag_value] = "D"

            if (L2_SIZE!=0 and L2_ASSOC!=0):
                L2_calling(L2_set_value,L2_tag_value,'r')



class sim:
    def __init__(self,BLOCKSIZE,L1_SIZE,L1_ASSOC,L2_SIZE,L2_ASSOC,REPLACEMENT_POLICY,INCLUSION_PROPERTY,trace_file):
        self.BLOCKSIZE=BLOCKSIZE
        self.L1_SIZE=L1_SIZE
        self.L1_ASSOC=L1_ASSOC
        self.L2_SIZE=L2_SIZE
        self.REPLACEMENT_POLICY=REPLACEMENT_POLICY
        self.INCLUSION_PROPERTY=INCLUSION_PROPERTY
        self.trace_file=trace_file

    def out(self):
        print("\n===== Simulator configuration =====")
        print('BLOCKSIZE:             ' + str(BLOCKSIZE))
        print('L1_SIZE:               '+str(L1_SIZE))
        print('L1_ASSOC:              '+str(L1_ASSOC))
        print('L2_SIZE:               '+str(L2_SIZE))
        print('L2_ASSOC:              '+str(L2_ASSOC))
        if REPLACEMENT_POLICY==0:
            print('REPLACEMENT POLICY:    LRU')
        if REPLACEMENT_POLICY==1:
            print('REPLACEMENT POLICY:    FIFO')
        if REPLACEMENT_POLICY==2:
            print('REPLACEMENT POLICY:    optimal')
        if INCLUSION_PROPERTY==1:
            print('INCLUSION PROPERTY:    inclusive')
        else:
            print('INCLUSION PROPERTY:    non-inclusive')
        print('trace_file:            '+str(trace_file))

        
obj=sim(BLOCKSIZE,L1_SIZE,L1_ASSOC,L2_SIZE,L2_ASSOC,REPLACEMENT_POLICY,INCLUSION_PROPERTY,trace_file)
obj.out()



print('===== L1 contents =====')
for i in range(len(L1)):
    if i<=9:
        text="Set     "+str(i)+":      "
    elif i<=99:
        text="Set     "+str(i)+":     "
    else:
        text="Set     "+str(i)+":    "
    for value in L1[i]:
        text += hex(int(value,2))[2:]
        if L1_dirty_dict[i][value]=='D':
            text+=' D  '
        else:
            text+='    '
    print(text)
if (L2_SIZE!=0 and L2_ASSOC!=0):
    print('===== L2 contents =====')
    for i in range(len(L2)):
        if i<=9:
            text="Set     "+str(i)+":      "
        elif i<=99:
            text="Set     "+str(i)+":     "
        else:
            text="Set     "+str(i)+":    "
        for value in L2[i]:
            text += hex(int(value,2))[2:]
            if L2_dirty_dict[i][value]=='D':
                text+=' D  '
            else:
                text+='    '
        print(text)
L1_miss_rate = float((L1_parameters['L1_readmiss']+L1_parameters['L1_writemiss'])/(L1_parameters['L1_reads']+L1_parameters['L1_writes']))
if (L2_SIZE!=0 and L2_ASSOC!=0):
    L2_miss_rate = float((L2_parameters['L2_readmiss'])/(L2_parameters['L2_reads']))
Traffic_number=0
print("===== Simulation results (raw) =====")
print("a. number of L1 reads:        "+str(L1_parameters['L1_reads']))
print("b. number of L1 read misses:  "+str(L1_parameters['L1_readmiss']))
print("c. number of L1 writes:       "+str(L1_parameters['L1_writes']))
print("d. number of L1 write misses: "+str(L1_parameters['L1_writemiss']))
print("e. L1 miss rate:              {:.6f}".format(L1_miss_rate))
print("f. number of L1 writebacks:   "+str(L1_parameters['L1_writebacks']))
print("g. number of L2 reads:        "+str(L2_parameters['L2_reads']))
print("h. number of L2 read misses:  "+str(L2_parameters['L2_readmiss']))
print("i. number of L2 writes:       "+str(L2_parameters['L2_writes']))
print("j. number of L2 write misses: "+str(L2_parameters['L2_writemiss']))
if (L2_SIZE!=0 and L2_ASSOC!=0):
    print("k. L2 miss rate:              {:.6f}".format(L2_miss_rate))
else:
    print("k. L2 miss rate:              "+str(0))
print("l. number of L2 writebacks:   "+str(L2_parameters['L2_writebacks']))
if L2_parameters['L2_reads']==0:
    Traffic_number = L1_parameters['L1_readmiss']+L1_parameters['L1_writemiss']+L1_parameters['L1_writebacks']
else:
    Traffic_number = L2_parameters['L2_readmiss']+L2_parameters['L2_writemiss']+L2_parameters['L2_writebacks']+inclusive_bit_value[0]
print("m. total memory traffic:      "+str(Traffic_number))