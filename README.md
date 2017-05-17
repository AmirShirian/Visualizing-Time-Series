# Online Spectral Image Generating for evolving cyclic Time Series
Visualize Your Data

Time series is one of the most important topic in scientific and financial application. A time series is a series of data in time order. The features of time series data include: large in size and high dimensions. Increasing usage of time series in science and industry, force scientist to attempt in this field. In this paper we demonstrate a new approach to visualize times series data based on frequency spectral. Moreover, it supports both cyclic and non-cyclic time series. This approach convert every time series to one picture that has its feature. This instrumentation can be used to reduced size, complexity and be raw data for a learning machine.
In this code I used AR model for the time serie and RLS algorithm to identify number of poles and zeros of a time serie in each time step. For more information you can read below paper.
http://rdcu.be/qY7m

For example when your data consist of some frequencies like below signal:

#S(t)=40cos(9(t−16))+50cos(13(t−26))−10cos(25(t−36))+80sin(3(t−46))+10

The output of this code may be:


![picture2](https://cloud.githubusercontent.com/assets/27130785/26144348/b25b5d6e-3afd-11e7-9e18-7aa1debd3e3b.png)


You can also see the variation of forgetting factor and number of modes:


![picture1](https://cloud.githubusercontent.com/assets/27130785/26144638/daa36a9a-3afe-11e7-8c8c-f2875b2d0dea.png)
