#%%
import pickle

# test_data = [['CACNA1S_matcher_combineFucn_testing1',["Unknown/Unknown"]],
#                 ['CACNA1S_matcher_combineFucn_testing2',["Reference/Reference"]]]#%%

test_data = {'CACNA1S_matcher_combineFucn_testing1' : ["Uown/Unknown"],
            'CACNA1S_matcher_combineFucn_testing2' : ["Refece/Reference"]}

def data_dumper(data=test_data):
    with open('expect_result.pickle', 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        
    # %%
if __name__ == '__main__':
    data_dumper()