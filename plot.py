# -*- coding: utf-8 -*-
"""
"""
# D:/Aghababaei/sampel/class_1/s289.csv



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# خواندن داده‌ها از فایل CSV
#file_path = 'H:/New folder/peakmi/y1000.csv'
file_path ='G:/Arshad/thesis/Aghili/beats/peakcontrol/s15545.csv'
data = pd.read_csv(file_path, header=None)

# بررسی داده‌ها
print(data.head())

# فرض می‌کنیم زمان نمونه‌برداری هر مقدار 1 میلی‌ثانیه است
sampling_rate = 1000  # نمونه در ثانیه (1 kHz)
num_samples = data.shape[0]
time = np.arange(num_samples) / sampling_rate  # بردار زمان

# استخراج سیگنال ECG (در اینجا تنها ستون اول را می‌گیریم)
ecg_signal = data[0]

# رسم سیگنال ECG
plt.figure(figsize=(10, 6))
plt.plot(time, ecg_signal, label='ECG Signal')
plt.xlabel('Time (s)')
