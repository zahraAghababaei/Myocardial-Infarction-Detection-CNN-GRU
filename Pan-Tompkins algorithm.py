# -*- coding: utf-8 -*-
"""
"""
import os
import numpy as np
import pandas as pd
import scipy.signal as signal

SN = 1
EN = 300

for k in range(1, 149):
    # مسیرها
    input_dir = "H:/New folder/signal  (" + str(k) + ")"
    output_dir = "H:/New folder/O1"
    
    # بررسی وجود مسیر خروجی و ایجاد آن در صورت عدم وجود
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # تابع Pan-Tompkins برای پردازش سیگنال ECG
    def pan_tompkins(ecg_signal, fs):
        # 1. Low pass filter
        b_low, a_low = signal.butter(1, 11 / (fs / 2), btype='low')
        ecg_low = signal.filtfilt(b_low, a_low, ecg_signal)
        
        # 2. High pass filter
        b_high, a_high = signal.butter(1, 5 / (fs / 2), btype='high')
        ecg_high = signal.filtfilt(b_high, a_high, ecg_low)
        
        # 3. Derivative filter
        derivative_filter = np.array([1, 2, 0, -2, -1]) * (fs / 8.0)
        ecg_derivative = np.convolve(ecg_high, derivative_filter, mode='same')
        
        # 4. Squaring
        ecg_squared = ecg_derivative ** 2
        
        # 5. Moving window integration
        window_size = int(0.150 * fs)
        integration_filter = np.ones(window_size) / window_size
        ecg_integrated = np.convolve(ecg_squared, integration_filter, mode='same')
        
        # 6. Thresholding
        threshold = 0.6 * np.max(ecg_integrated)
        detected_peaks = np.where(ecg_integrated > threshold)[0]
        
        # شناسایی پیک‌های QRS
        peaks, _ = signal.find_peaks(ecg_integrated, distance=fs/2, height=threshold)
        
        return peaks
    
    # دریافت لیست فایل‌های CSV در مسیر ورودی
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    
    # شمارنده برای نام‌گذاری نمونه‌ها
    sample_counter = SN
    
    # پردازش هر فایل CSV
    for csv_file in csv_files:
        file_path = os.path.join(input_dir, csv_file)
        
        # خواندن داده‌های CSV
        data = pd.read_csv(file_path)
        
        # فرض می‌کنیم داده‌های ECG در ستون اول هستند
        ecg_signal = data.iloc[:, 0].values
        
        # تعیین نرخ نمونه‌برداری
        sampling_rate = 250  # Hz
        
        # استفاده از الگوریتم Pan-Tompkins برای تشخیص پیک‌های QRS
        r_peaks = pan_tompkins(ecg_signal, sampling_rate)
        
        # صرف نظر کردن از دو تناوب اول
        if len(r_peaks) > 2:
            r_peaks = r_peaks[2:]
            r_peaks = r_peaks[~np.isnan(r_peaks)]
        
        # محاسبه فاصله R-R
        rr_intervals = np.diff(r_peaks)
        rr_intervals = rr_intervals[~np.isnan(rr_intervals)]
        avg_rr_interval = np.mean(rr_intervals) if len(rr_intervals) > 0 else 0
        
        # استخراج قطعات PQRST بر اساس نقاط R
        for r_peak in r_peaks:
            if np.isnan(r_peak) or avg_rr_interval == 0:
                continue
            
            start = int(r_peak - 0.2 * avg_rr_interval)
            end = int(r_peak + 0.4 * avg_rr_interval)
            
            # اطمینان از اینکه محدوده در محدوده سیگنال باشد
            start = max(start, 0)
            end = min(end, len(ecg_signal))
            
            if start >= end:
                continue
            
            pqrst_segment = ecg_signal[start:end]
            
            # تبدیل مقدار نمونه به DataFrame
            sample_df = pd.DataFrame(pqrst_segment)
            sample_file = os.path.join(output_dir, f'y{sample_counter}.csv')
            sample_df.to_csv(sample_file, index=False, header=False)
            sample_counter += 1
            
            # توقف پس از ذخیره 300 نمونه
            if sample_counter > EN:
                break
        
        # توقف حلقه اصلی در صورت ذخیره 300 نمونه
        if sample_counter > EN:
            break
    
    print("k=" + str(k))
    SN += 300
    EN += 300
