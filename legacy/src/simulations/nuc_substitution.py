'''
Created on Oct 30, 2013

@author: marieidowu
'''
from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
import matplotlib.figure as fig
import numpy as np
import random as rand
from math import fabs
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def get_G(A,C,T):
    ''' 
    returns the value of G ensuring A,C,T,G sum to 1
    '''
        
    G = fabs(1-(A+C+T))
    
    return G
def validate_A(A,C,T):
    '''
    Checks that the value of A is greater than or equql to zero
    and less than 1-C-T
    returns the value of A and G ensuring A,C,T,G sum to 1
    '''
    a_min =0.0
    a_max = 1-C-T
    
    if A < a_min:
        A = a_min
    elif A > a_max:
        A = a_max
        
    G = fabs(1-(A+C+T))
    
    return A,G
        
def validate_T(A,C,T):
    '''
    Checks that the value of T is greater than or equql to zero
    and less than 1-C-A
    returns the value of T and G ensuring A,C,T,G sum to 1
    '''
    a_min =0.0
    a_max = 1-C-A
    
    if T < a_min:
        T = a_min
    elif T > a_max:
        T = a_max
        
    G = fabs(1-(A+C+T))
    
    return T,G

def validate_C(A,C,T):
    '''
    Checks that the value of C is greater than or equql to zero
    and less than 1-T-A
    returns the value of C and G ensuring A,C,T,G sum to 1
    '''
    a_min =0.0
    a_max = 1-T-A
    
    if C < a_min:
        C = a_min
    elif C > a_max:
        C = a_max
         
        
    G = fabs(1-(A+C+T))
    
    return C,G

def validate_transition(value):
    '''
    Set rate of transitions, these are assumed to be equal under time-reversible model
    Transitions are:
        C->T and T->C
        A->C and C->A
        A->T and T->A
        G->C and C->G
        G->T and T->G
        
    '''
    
    min_transition = 0.000001
    max_transition = 0.5
    
    if value < min_transition:
        value = min_transition
    elif value > max_transition:
        value = max_transition
     
    return value

def base_frequencies(sequence_1, sequence_2, sequence_length):
    '''
    get observed base frequencies for a sequence. Both base pairs must be
    identical to be counted. See table 3.1 in Nei and Kumar 2000. 
    DNA sequences are coded as A=1 T=2 C=3 G=4
    '''
    observed_n_A, observed_n_T,observed_n_C,observed_n_G,n = 0,0,0,0,0
    
    for base_pair in xrange(sequence_length):
        sq_1_bp= sequence_1[0, base_pair]
        sq_2_bp= sequence_2[0, base_pair]
       
        if sq_1_bp == 1 and sq_2_bp==1:
            observed_n_A += 1
            n+=1
            #print 'n should increase by 1 %d' %n
        elif sq_1_bp == 2 and sq_2_bp==2:
            observed_n_T += 1
            n+=1
            #print 'n should increase by 1 %d' %n
        elif sq_1_bp == 3 and sq_2_bp==3:
            observed_n_C += 1
            n+=1
            #print 'n should increase by 1 %d' %n
        elif sq_1_bp == 4 and sq_2_bp==4:
            observed_n_G += 1
            n+=1
            #print 'n should increase by 1 %d' %n
    #print 'what is n? %d' %n
    p_A = observed_n_A/n
    p_T = observed_n_T/n
    p_C = observed_n_C/n
    p_G = observed_n_G/n
    return p_A, p_T, p_C, p_G

def compare_sequences(sequence_1,sequence_2,sequence_length):
    '''
    compare sequences to count up differences of all types.
    DNA sequences are coded as A=1 T=2 C=3 G=4
    '''
    AG,CT,AC,AT,GC,GT = 0,0,0,0,0,0
    
    for base_pair in xrange(sequence_length):
        seq1_entry = sequence_1[0, base_pair]
        seq2_entry = sequence_2[0,base_pair]
        if seq1_entry == 1 and seq2_entry==2:
            AT = AT+1
        elif seq1_entry == 1 and seq2_entry==3:
            AC = AC+1
        elif seq1_entry == 1 and seq2_entry==4:
            AG = AG+1
        elif seq1_entry == 2 and seq2_entry==1:
            AT = AT+1
        elif seq1_entry == 2 and seq2_entry==3:
            CT = CT+1
        elif seq1_entry == 2 and seq2_entry==4:
            GT = GT+1
        elif seq1_entry == 3 and seq2_entry==1:
            AC = AC+1
        elif seq1_entry == 3 and seq2_entry==2:
            CT = CT+1
        elif seq1_entry == 3 and seq2_entry==4:
            GC = GC+1
        elif seq1_entry == 4 and seq2_entry==1:
            AG = AG+1
        elif seq1_entry == 4 and seq2_entry==2:
            GT = GT+1
        elif seq1_entry == 4 and seq2_entry==3:
            GC = GC+1
            
    return AG,CT,AC,AT,GC,GT
        
def mutate_sequence(sequence,base_pair,AG,CT,AC,AT,GC,GT):
    '''
    determine state change at a given base pair j in sequence "seq" when there has been a
    mutation. Also, count up the actual mutational changes that occur. 
    DNA sequences are coded as A=1 T=2 C=3 G=4.
    '''
    
    if sequence[0,base_pair] == 1: #base is A
        probability = rand.random()*(AG+AC+AT)
        if probability <= AG: #mutation is A->G transition
            sequence[0,base_pair]=4
            type_mutation = 1
        elif probability > AG and probability < AC+AG: #mutation is A->C transversion
            sequence[0,base_pair]=3
            type_mutation = 3
        else: #mutation is A->T transversion
            sequence[0,base_pair] = 2
            type_mutation=4
    
    elif sequence[0,base_pair]==2: #base is T
        probability = rand.random()*(CT+AT+GT)
        if probability <= CT: #mutation is T->C transition
            sequence[0,base_pair]=3
            type_mutation = 2
        elif probability > CT and probability < AT+CT: #mutation is T->A transversion
            sequence[0,base_pair]=1
            type_mutation = 4
        else: #mutation is T->G transversion
            sequence[0,base_pair] = 4
            type_mutation=6
   
    elif sequence[0,base_pair]==3: #base is C
        probability = rand.random()*(CT+AC+GC)
        if probability <= CT: #mutation is C->T transition
            sequence[0,base_pair]=2
            type_mutation = 2
        elif probability > CT and probability < AC+CT: #mutation is C->A transversion
            sequence[0,base_pair]=1
            type_mutation = 3
        else: #mutation is C->G transversion
            sequence[0,base_pair] = 4
            type_mutation=5
    
    elif sequence[0,base_pair]==4: #base is G
        probability = rand.random()*(AG+GC+GT)
        if probability <= AG: #mutation is G->A transition
            sequence[0,base_pair]=1
            type_mutation = 5
        elif probability > AG and probability < GC+AG: #mutation is G->C transversion
            sequence[0,base_pair]=3
            type_mutation = 5
        else: #mutation is G->T transversion
            sequence[0,base_pair] = 4
            type_mutation=6
            
    return sequence, type_mutation

def prepare_plot_data(A,T,C,G,AG,CT,AC,AT,GC,GT,sequence_length,generations):
    
    sequence_1 = create_DNA_sequence(sequence_length,A,T,C)
    sequence_2 = sequence_1.copy()
    probability_mutation = AG+CT+AT+GC+GT
    true_AC,true_AG,true_AT,true_CT,true_GC,true_GT,nucleotide_AG,nucleotide_CT, nucleotide_AC,nucleotide_AT,nucleotide_GC,nucleotide_GT, p_A, p_T,p_C,p_G =get_true_nucleotide_differences(generations,sequence_length,probability_mutation, sequence_1,sequence_2,AG,CT,AC,AT,GC,GT)
    
    p_observed_differences = (nucleotide_AG+nucleotide_CT+nucleotide_AC+nucleotide_AT+nucleotide_GC+nucleotide_GT)/sequence_length # element wise division
    P = (nucleotide_AG+nucleotide_CT)/sequence_length
    Q = (nucleotide_AC+nucleotide_AT+nucleotide_GC+nucleotide_GT)/sequence_length
    
    
    total_true_differences =  true_AG + true_CT + true_AC + true_AT + true_GC + true_GT
    total_true_transitions = true_AG + true_CT
    total_true_transversions = true_AC + true_AT + true_GC + true_GT
    
    cumulative_true_differences = np.cumsum(total_true_differences, axis=1)
#    #express total true changes, total true transversions, and total true transitions 
    #per nucleotide site for each sequence
    p_cumulative_true_differences = cumulative_true_differences/sequence_length
#    #sum per site total true changes for the two sequences (sum is across columns)
    p_cumulative_true_differences_both = np.sum(p_cumulative_true_differences, axis=0)
   
    #get distance estimate using Jukes-Cantor model
    jukes_cantor_distance = -0.75*(np.log(1-1.25*p_observed_differences))

    kimura_1980_distance = -0.5*np.log(1-2*P-Q) - 0.25*np.log(1-2*Q)
    #get distance estimate using Tamura 1992 model
    CG_content = (p_C+p_G)
    h = 2*CG_content*(1-CG_content)

    tamura_1992_distance = -h*np.log(1-(P/h)- Q) - 0.5 *(1-h)*np.log(1-2*Q)

    #get distance estimate using Tamura & Nei 1993 model
    
    p_R = p_A + p_G    # observed freq of purines
    p_Y = p_T + p_C    # observed freq of pyrimidines
    P1 = nucleotide_AG/sequence_length  # A-G transitions
    P2 = nucleotide_CT/sequence_length # T-C transitions
    
    # see http://www.megasoftware.net/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/nucleotide_substitution_models/hc_tamura_nei_distance.htm
    k1 = (2*p_A*p_G)/p_R
    k2 = (2*p_T*p_C)/p_Y
    k3 = 2*(p_R*p_Y - ((p_A*p_G*p_Y)/p_R) - (p_T*p_C*p_R)/p_Y)
#    k4 = 2*(p_A*p_G + p_T*p_C + p_R*p_Y)
    
    w1 = 1 - P1/k1 - Q/(2*p_R)
    w2 = 1 - P2/k2 - Q/(2*p_Y)
    w3 = 1 - Q/(2*p_R*p_Y)
#    
#    print 'w1 %s' %w1
#    print 'w2 %s' %w2
#    print 'w3 %s' %w3
    w1,w2,w3, generations_until_negative = negative_value_test(w1, w2, w3, generations)
    tamura_nei_1993_distance = -k1*np.log(w1) - k2*np.log(w2) - k3*np.log(w3)
    
    return p_observed_differences, p_cumulative_true_differences_both, jukes_cantor_distance,kimura_1980_distance,tamura_1992_distance,tamura_nei_1993_distance,generations_until_negative

def plot(A,T,C,G,AG,CT,AC,AT,GC,GT,sequence_length,generations, clear):
    figure = fig.Figure((16,13),facecolor='#ececec')
    axes = figure.add_subplot(1,1,1)
    figure.suptitle('Apparent and nucleotide substitution model corrected DNA sequence divergence with time', size=18,color='#003E7E',family='serif', weight='bold')

    if not clear:
    	p_observed_differences, p_cumulative_true_differences_both, jukes_cantor_distance,kimura_1980_distance,tamura_1992_distance,tamura_nei_1993_distance, generations_until_negative = prepare_plot_data(A, T, C, G, AG, CT, AC, AT, GC, GT, sequence_length, generations)

       	t92_generations_until_unity, tn93_range = check_values_greater_than_one(tamura_nei_1993_distance, tamura_1992_distance, generations, generations_until_negative)
        axes.set_xlim((1,generations+1))
        axes.set_ylim((0.0,1.05))
        axes.set_ylabel("Sequence divergence (proportion of sites)", fontdict={'size':18,'color':'#002147','family':'serif'})
        axes.set_xlabel("Generation", fontdict={'size':18,'color':'#002147','family':'serif'})
        
        axes.plot(p_observed_differences[0,:],label='Observed divergence') #labels are used for the legend
        axes.plot(p_cumulative_true_differences_both,linestyle='--',color='k', label='Sites experiencing mutation in both sequences')
        axes.plot(jukes_cantor_distance[0,:],linestyle='-.',color='r', label='Juke-Cantor estimate')         
        axes.plot(kimura_1980_distance[0,:],color='g', label='Kimura 1980 estimate')
        axes.plot(tamura_1992_distance[0,0:t92_generations_until_unity],color='m', label='Tamura 1992 estimate')
        axes.plot(tamura_nei_1993_distance[0,0:tn93_range],color='k', label='Tamura & Nei 1993 estimate')
        axes.legend(loc='best',prop={'size':18})
    canvas = FigureCanvas(figure)
    return canvas

def negative_value_test(w1,w2,w3,generations):
    '''
    test w1, w2, and w3 for negative values that will give imaginary result if log(w)
    assign limit_value in those cases to prevent taking a log of zero or a negative. 
    '''
    
    generations_until_negative = generations
    limit_value = 0.01
    
    for i in xrange( generations-1,-1,-1): #xrange([start,] stop[, step])
        if w1[0,i] < 0:
            w1[0,i] = limit_value
            generations_until_negative = i
        
        if w2[0,i] < 0:
            w2[0,i] = limit_value
            generations_until_negative = i
        
        if w3[0,i] < 0:
            w3[0,i] = limit_value
            generations_until_negative = i
    return w1,w2,w3, generations_until_negative

def check_values_greater_than_one(tn93_distance, t92_distance, generations,generations_until_negative):
    '''
    check for values of T92_dist and TN93_dist greater than 1.0
    '''
    generations_until_unity = generations
    t92_generations_until_unity = generations
    
    for i in xrange(generations-1,-1,-1):
        if tn93_distance[0,i] > 1.0:
            tn93_generations_until_unity =i
        
        if t92_distance[0,i] >1.0:
            t92_generations_until_unity = i
    
    if generations_until_unity > generations_until_negative:
        tn93_range = generations_until_negative
    else:
        tn93_range = generations_until_unity
                
    return t92_generations_until_unity, tn93_range

def init_base_pair_frequencies(generations):

    '''
    allocate vectors to store 
    '''
    vector = np.zeros((1,generations))
    return vector, vector.copy(), vector.copy(), vector.copy()
        
def get_true_nucleotide_differences(generations,sequence_length,probability_mutation, sequence_1,sequence_2,AG,CT,AC,AT,GC,GT):
    
    '''
    allocate vectors to store true nucleotide differences at each
    generation. First row is for sequence_1 second row is for sequence_2, columns are
    for generations
    '''
    nucleotide_AG,nucleotide_CT, nucleotide_AC,nucleotide_AT,nucleotide_GC,nucleotide_GT = init_nucleotide_vectors(generations,1)
    true_AG,true_CT, true_AC,true_AT,true_GC,true_GT = init_nucleotide_vectors(generations, 2)
    p_A, p_T,p_C,p_G = init_base_pair_frequencies(generations)
    
    for i in xrange(generations):
        for j in xrange(sequence_length):
            probability = rand.random()
            if probability > 1- probability_mutation:
                sequence_1, cat_change = mutate_sequence(sequence_1, j, AG, CT, AC, AT, GC, GT)
                true_AC,true_AG,true_AT,true_CT,true_GC,true_GT = populate_true_nucleotide_differences(cat_change,true_AC,true_AG,true_AT,true_CT,true_GC,true_GT,0,i)
            probability = rand.random()
            if probability > 1- probability_mutation:
                sequence_2, cat_change = mutate_sequence(sequence_2, j, AG, CT, AC, AT, GC, GT)
                true_AC,true_AG,true_AT,true_CT,true_GC,true_GT = populate_true_nucleotide_differences(cat_change,true_AC,true_AG,true_AT,true_CT,true_GC,true_GT,1,i)
        nucleotide_AG[0,i],nucleotide_CT[0,i], nucleotide_AC[0,i],nucleotide_AT[0,i],nucleotide_GC[0,i],nucleotide_GT[0,i] = compare_sequences(sequence_1, sequence_2, sequence_length)
        p_A[0,i],p_T[0,i],p_C[0,i],p_G[0,i] = base_frequencies(sequence_1, sequence_2, sequence_length)
        
    return true_AC,true_AG,true_AT,true_CT,true_GC,true_GT,nucleotide_AG,nucleotide_CT, nucleotide_AC,nucleotide_AT,nucleotide_GC,nucleotide_GT, p_A, p_T,p_C,p_G

def populate_true_nucleotide_differences(cat_change,true_AC,true_AG,true_AT,true_CT,true_GC,true_GT, row,column):
    if cat_change ==1:
        true_AG[row,column] = true_AG[row,column]+1
    elif cat_change ==2:
        true_CT[row,column] = true_CT[row,column]+1
    elif cat_change ==3:
        true_AC[row,column] = true_AC[row,column]+1
    elif cat_change ==4:
        true_AT[row,column] = true_AT[row,column]+1
    elif cat_change ==5:
        true_GC[row,column] = true_GC[row,column]+1
    else:
        true_GT[row,column] = true_GT[row,column]+1
    
    return true_AC,true_AG,true_AT,true_CT,true_GC,true_GT
        
def init_nucleotide_vectors(generations,rows):
    '''
    allocate vectors to store observed nucleotide differences at each
    generation
    '''
    vector = np.zeros((rows, generations))
    return vector, vector.copy(), vector.copy(), vector.copy(), vector.copy(), vector.copy()
    
def create_DNA_sequence(sequence_length, A,T,C):
    sequence_1 = np.zeros((1,sequence_length))
    
    for i in xrange(sequence_length):
        probability = rand.random()
        
        if probability <= A:
            sequence_1[0,i]=1
        elif probability > A and probability <= A+T:
            sequence_1[0,i]=2
        elif probability > A+T and probability <= A+T+C:
            sequence_1[0,i]=3
        else: #only other option is G
            sequence_1[0,i]=4
        
    return sequence_1
            