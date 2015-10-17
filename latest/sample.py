import numpy as np
import pylab as plt
m=1000
p=[]
st=[]
s_f =float(1e-12)
n=4096


for std in range(1,3):
    noise_fft=np.zeros(n/2)
    st.append(std)
    for i in xrange(m):
        noise=np.random.standard_normal(n)*std
        noise_fft += np.abs(np.fft.fft(noise,n)[0:n/2])**2
    plt.plot(noise);plt.show()
    plt.plot(noise_fft);plt.show()
    q=np.mean(noise_fft)
    print std,"\t",q
    p.append(q)
plt.plot(st,p);plt.show()
