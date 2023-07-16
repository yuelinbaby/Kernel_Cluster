from keras.layers import Input, Dense
from keras.models import Model
from keras.datasets import mnist
import numpy as np
import pandas as pd
import sys

import argparse
parser = argparse.ArgumentParser(description='para transfer')
parser.add_argument('--indir', type=str, help='indir -> inputdir.')
parser.add_argument('--samid', type=str, help='samid -> sample id.')
#parser.add_argument('--downsam', type=str, help='downsam -> downsampleratio')
args = parser.parse_args()

indir="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output"
samid="151673"

samdir = args.indir
samid = args.samid

def encode_kernel_pcf(samdir,k,p):
    #tagfile = f'{samdir}/{samid}.ground_truth.csv'
    expfile = f'{samdir}/kernel_{k}_{p}.exp_data.csv'

    # data preprocessing
    #tag = pd.read_csv(tagfile,index_col=0)
    exp_ori = pd.read_csv(expfile,index_col=0)

    #MinMax Scaling

    from sklearn.preprocessing import MinMaxScaler

  
    scaler = MinMaxScaler() #axis=1 

    exp = scaler.fit_transform(exp_ori.T)
    exp=exp.T


    # 
    samnum=exp.shape[1]
    genenum=exp.shape[0]
    arr = np.arange(samnum)

    # 
    trainnum=int(samnum*0.8)
    sample_train = np.random.choice(arr, size=trainnum, replace=False)
    sample_train = np.sort (sample_train)

    sample_test = np.setdiff1d(arr,sample_train) 
    sample_test = np.sort(sample_test)

    # datapreprocessing
    exp2 = exp
    x_train = exp2[:,sample_train].T # (3380, 15412)
    x_test = exp2[:,sample_test].T # (846, 15412)

    x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
    x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))

    # 定义自编码器的结构
    input_img = Input(shape=(genenum,))
    encoded = Dense(32, activation='relu')(input_img)
    decoded = Dense(genenum, activation='sigmoid')(encoded)

    # 建立自编码器模型
    autoencoder = Model(input_img, decoded)

    # 编译模型
    autoencoder.compile(optimizer='adam', loss='mse')

    # 训练模型
    autoencoder.fit(x_train, x_train,
                    epochs=100,
                    batch_size=256,
                    shuffle=True,
                    validation_data=(x_test, x_test))

    # 评估模型
    mse = autoencoder.evaluate(x_test, x_test)
    print('MSE: %.4f' % mse)
    return(mse)


#ks = [3,5,7,9,11]
ks = [3]
ps = list(range(0,9))

mse_m=np.zeros((len(ks),len(ps)))
mseoutfile = f'{samdir}/{samid}.encorder.mse.csv'
for k in ks:
    for p in ps:
        mset = encode_kernel_pcf(samdir,k,p)
        k2=int((k-1)/2-1)
        p2=int(p-1)
        mse_m[k2][p]=mset
np.savetxt(mseoutfile,mse_m)

    



