library(dplyr)
library(zoo)  # 이동평균 계산
library(ggplot2)
library(reshape2)
library(lubridate)
library(tidyr)
library(moments)  # 첨도 계산을 위한 패키지
library(forecast)  # ACF 분석
library(TSA)       # 주파수 분석(FFT)
library(signal)  # 신호 필터링
library(gridExtra)  # 여러 개의 ggplot을 한 화면에 표시
library(grid)

load("data//run1.RData")
load("data//run2.RData")
load("data//run3.RData")
load("data//run4.RData")


# 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)
yearSet = year(dates) %>% unique()
years = year(dates)
# 연도별 0~1 변환
df_time <- data.frame(
  Date = dates,
  Year = year(dates),
  DayOfYear = yday(dates),  # 연도의 몇 번째 날인지
  TotalDays = ifelse(leap_year(dates), 366, 365)  # 윤년 고려하여 총 일수 설정
)

# 변환된 시간 계산 (연도를 소수점으로 변환)
df_time <- df_time %>%
  mutate(Time_Normalized = (Year-1850) + (DayOfYear / TotalDays))




# 데이터 선택 (1행 3열의 시계열 데이터)
ts_data <- run1[2, 3, ]

# FFT 변환 수행
fft_result <- fft(ts_data)
power_spectrum <- Mod(fft_result)^2  # 진폭의 제곱 (강도)
n <- length(ts_data)  # 데이터 길이

# 주파수 및 주기 계산
freqs <- (0:(n-1)) / n  # 주파수 값
periods <- 1 / freqs    # 주기 (일 단위)
valid_idx <- 2:(n/2)    # DC 성분(0번 주파수) 제거, 대칭이라 절반만 사용

# 데이터 프레임 생성
fft_df <- data.frame(Period = periods[valid_idx], Power = power_spectrum[valid_idx])

# 로그 스케일을 적용하여 그래프 그리기
ggplot(fft_df, aes(x = Period, y = Power)) +
  geom_point(color = "blue", size = 2) +                     # 점 그래프
  geom_line(linetype = "dashed", color = "blue", alpha = 0.7) +  # 선 연결
  scale_x_log10() +  # 로그 스케일 (주기 간 차이 큼)
  labs(title = "FFT 분석 결과 - 주기별 강도",
       x = "주기 (일)", 
       y = "강도") +
  theme_minimal()

365/4
fft_df %>% filter(Power>210000)








# 예제 데이터 (run1이 주어진 상태라고 가정)
# run1은 (5,5,60225) 형태의 다차원 배열

# 1. 5x5 지역의 합 벡터 계산
total_sum <- apply(run1, 3, sum)  # (60225,) 형태

# 2. (1,1) 위치의 강수량 데이터
target_series <- run1[1,1,]  # (60225,) 형태

# 3. 상관계수 계산
correlation <- cor(total_sum, target_series, method = "pearson")

# 결과 출력
print(correlation)





# 특정 위치의 시계열 데이터 선택
ts_data <- run1[5, 2, ]

# 총 기간 (60225일)과 연도 개수 계산
num_years <- floor(length(ts_data) / 365)  # 전체 연도 수

# 매년 두 번째 날짜(1월 2일) 데이터 추출
yearly_data <- ts_data[seq(10, num_years * 365, by = 365)]

ccf(yearly_data, yearly_data, lag.max = 10, main = "Cross-Correlation (X_t vs X_{t-lag})")







library(ggplot2)
library(reshape2)

# 5x5 지역의 시간 데이터 추출
data_matrix <- matrix(NA, nrow = 60225, ncol = 25)

# 5x5 지역 데이터를 열 단위로 저장
index <- 1
for (i in 1:5) {
  for (j in 1:5) {
    data_matrix[, index] <- run1[i, j, ]
    index <- index + 1
  }
}

# 열 이름 설정 (위치 좌표)
colnames(data_matrix) <- paste0("(", rep(1:5, each = 5), ",", rep(1:5, times = 5), ")")

# 상관행렬 계산
cor_matrix <- cor(data_matrix, use = "complete.obs")

# 데이터 프레임 변환 (히트맵을 위해 long format으로 변환)
melted_cor <- melt(cor_matrix)

# 히트맵 시각화
ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
  labs(title = "5x5 지역 간 상관계수 (Correlation Matrix)", x = "지역 좌표", y = "지역 좌표", fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 5x5 지역의 시간 데이터 추출
data_matrix <- matrix(NA, nrow = 60225, ncol = 25)

# 5x5 지역 데이터를 열 단위로 저장
index <- 1
for (i in 1:5) {
  for (j in 1:5) {
    data_matrix[, index] <- run1[i, j, ]
    index <- index + 1
  }
}

# 공분산 행렬 계산 (분산-공분산 행렬)
var_matrix <- var(data_matrix, use = "complete.obs")

# 5x5 형태로 변환 (대각선 값 = 각 지역의 분산)
var_5x5 <- matrix(diag(var_matrix), nrow = 5, ncol = 5)

# 결과 출력
print(var_5x5)











# 특정 위치의 시계열 데이터 선택
ts_data <- run1[2, 3, ]  # 특정 (1,3) 위치 선택

# 총 연도 개수 계산
num_years <- floor(length(ts_data) / 365)  # 전체 연도 수

# 1~365일별 분산을 저장할 벡터
variance_values <- numeric(365)

# 각 날짜(1~365일)에 대해 분산 계산
for (i in 1:365) {
  daily_values <- ts_data[seq(i, num_years * 365, by = 365)]  # 매년 같은 날짜 값 추출
  variance_values[i] <- var(daily_values, na.rm = TRUE)  # 분산 계산
}

# 데이터 프레임 생성
df <- data.frame(Day = 1:365, Variance = variance_values)

# 시각화
library(ggplot2)
ggplot(df, aes(x = Day, y = Variance)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1) +
  labs(title = "Variance of { X_{i+365*n} } by Day",
       x = "Day (1~365)",
       y = "Variance") +
  theme_minimal()




vv = 365+1:365
# 5x5 지역의 상관행렬을 저장할 행렬
cor_matrix <- matrix(NA, nrow = 5, ncol = 5)

# 기준 셀 (3,3) 위치의 시계열 데이터
reference_series <- run1[3,3,vv]

# 각 셀과 (3,3) 위치의 상관계수 계산
for (i in 1:5) {
  for (j in 1:5) {
    cor_matrix[i, j] <- cor(run1[i, j, vv], reference_series, use = "complete.obs")
  }
}

# 결과 출력
print(cor_matrix)




# 특정 위치의 시계열 데이터 선택
ts_data <- run1[1, 3, ]

# 총 연도 개수 계산
num_years <- floor(length(ts_data) / 365)

# 월별 분산을 저장할 벡터
variance_values <- numeric(12)

# x축 월 정보
months <- 1:12

# 각 월별 분산 계산
for (month in 1:12) {
  # 월별 해당 날짜(연도별) 추출 (예: 1월 = 1, 32, 63, ..., 12월 = 335, 366, ...)
  month_indices <- seq(month, num_years * 365, by = 365)
  monthly_values <- ts_data[month_indices]  # 월별 데이터
  variance_values[month] <- var(monthly_values, na.rm = TRUE)  # 분산 계산
}

# 데이터 프레임 생성
df <- data.frame(Month = months, Variance = variance_values)

# 시각화
library(ggplot2)
ggplot(df, aes(x = Month, y = Variance)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # x축을 1~12 (Jan, Feb, ...)로 표시
  labs(title = "Monthly Variance",
       x = "Month",
       y = "Variance") +
  theme_minimal()



# 필요한 라이브러리 로드
library(ggplot2)
library(dplyr)
library(lubridate)
library(reshape2)

# 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)

# 월만 추출 (1~12월)
months <- month(dates)

# 5x5 지역 데이터를 long-format 데이터 프레임으로 변환
data_list <- list()

for (i in 1:5) {
  for (j in 1:5) {
    location_name <- paste0("(", i, ",", j, ")")  # 위치명을 문자열로 생성
    df <- data.frame(Month = months, Value = run1[i, j, ], Location = location_name)
    data_list[[location_name]] <- df
  }
}

# 하나의 데이터 프레임으로 병합
df <- do.call(rbind, data_list)

# 월별 평균값 계산 (연도를 무시하고 1~12월 기준으로 그룹화)
df_monthly <- df %>%
  group_by(Month, Location) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# ggplot을 사용하여 5x5 지역별 월별 평균 시각화
ggplot(df_monthly, aes(x = Month, y = Value, color = Location)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # x축을 "Jan, Feb, ..." 형식으로 변환
  labs(title = "Monthly Average for 5x5 Locations",
       x = "Month",
       y = "Average Value",
       color = "Location") +
  theme_minimal()














# 특정 위치의 시계열 데이터 선택
ts_data <- run1[5,1, ]

# 총 연도 개수 계산
num_years <- floor(length(ts_data) / 365)

# 월별 분산을 저장할 벡터
variance_values <- numeric(12)

# x축 월 정보
months <- 1:12

# 각 월별 분산 계산
for (month in 1:12) {
  # 월별 해당 날짜(연도별) 추출 (예: 1월 = 1, 32, 63, ..., 12월 = 335, 366, ...)
  month_indices <- seq(month, num_years * 365, by = 365)
  monthly_values <- ts_data[month_indices]  # 월별 데이터
  variance_values[month] <- var(monthly_values, na.rm = TRUE)  # 분산 계산
}

# 데이터 프레임 생성
df <- data.frame(Month = months, Variance = variance_values)

# 시각화
library(ggplot2)
ggplot(df, aes(x = Month, y = Variance)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # x축을 1~12 (Jan, Feb, ...)로 표시
  labs(title = "Monthly Variance",
       x = "Month",
       y = "Variance") +
  theme_minimal()










# 필요한 라이브러리 로드
library(ggplot2)
library(dplyr)
library(lubridate)

# 특정 위치의 시계열 데이터 선택
ts_data <- run1[1, 3, ]

# 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)

# 연도 및 월 정보 추출
years <- year(dates)
months <- month(dates)

# 데이터 프레임 생성
df <- data.frame(Year = years, Month = months, Value = ts_data)

# 연도별 월별 분산 계산
df_variance <- df %>%
  group_by(Year, Month) %>%
  summarise(Variance = var(Value, na.rm = TRUE), .groups = "drop")

# ggplot을 사용하여 연도별 월별 분산 시각화
ggplot(df_variance[1:200,], aes(x = Month, y = Variance)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # x축을 "Jan, Feb, ..." 형식으로 변환
  facet_wrap(~ Year, scales = "free_y") +  # 연도별로 분할
  labs(title = "Monthly Variance by Year",
       x = "Month",
       y = "Variance") +
  theme_minimal()

                     






library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)

# 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)

# 모든 위치(5x5)에서 데이터를 추출하여 데이터 프레임 생성
location_list <- list()
for (i in 1:5) {
  for (j in 1:5) {
    loc_name <- paste0("Location_", i, "_", j)  # 위치 이름 지정
    location_list[[loc_name]] <- run1[i, j, ]  # 각 위치의 데이터 저장
  }
}

# 데이터 프레임 생성 (각 위치가 열로 구성됨)
df <- data.frame(Date = dates, location_list)

# 연도 및 월 정보 추가
df <- df %>%
  mutate(
    Year = year(Date),
    Month = month(Date)
  )

# 데이터 long format으로 변환 (tidy format)
df_long <- df %>%
  pivot_longer(cols = starts_with("Location"), 
               names_to = "Location", 
               values_to = "Value")

# 연도별 월별 
df_variance <- df_long %>%
  group_by(Year, Month, Location) %>%
  summarise(Median = quantile(Value,0.9, na.rm = TRUE), .groups = "drop")

# ggplot을 사용하여 연도별 월별 분산 시각화 (25개 위치)
ggplot(df_variance[1:4000, ], 
       aes(x = Month, y = Median, color = Location, group = Location)) +
  geom_line(size = 0.8) +   # 선 그래프
  geom_point(size = 1) +  # 데이터 포인트 추가
  scale_x_continuous(breaks = 1:12, labels = month.abb) +  # x축을 "Jan, Feb, ..." 형식으로 변환
  facet_wrap(~ Year, scales = "free_y", ncol = 5) +  # 연도별로 분할, 5열로 배치
  labs(title = "Monthly Variance by Year and Location (5×5 Grid)",
       x = "Month",
       y = "Variance",
       color = "Location") +
  theme_minimal()











library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)

# 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)

# 모든 위치(5x5)에서 데이터를 추출하여 데이터 프레임 생성
grid_data <- expand.grid(X = 1:5, Y = 1:5)  # 5x5 그리드 좌표 생성
grid_data$Location <- paste0("Location_", grid_data$X, "_", grid_data$Y)  # 위치 이름

# 각 위치의 데이터를 리스트로 저장
location_list <- list()
for (i in 1:5) {
  for (j in 1:5) {
    loc_name <- paste0("Location_", i, "_", j)  # 위치 이름 지정
    location_list[[loc_name]] <- run1[i, j, ]  # 각 위치의 데이터 저장
  }
}

# 데이터 프레임 생성 (각 위치가 열로 구성됨)
df <- data.frame(Date = dates, location_list)

# 연도 정보 추가
df <- df %>%
  mutate(Year = year(Date))

# 데이터 long format으로 변환 (tidy format)
df_long <- df %>%
  pivot_longer(cols = starts_with("Location"), 
               names_to = "Location", 
               values_to = "Value")

# 연도별 90% 분위수 계산
df_quantile90 <- df_long %>%
  group_by(Year, Location) %>%
  # summarise(Quantile_90 = quantile(Value, 0.9, na.rm = TRUE), .groups = "drop")
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")
# 5x5 grid 좌표를 df_quantile90과 병합
df_heatmap <- left_join(df_quantile90, grid_data, by = "Location")

# 연도별 5×5 그리드 히트맵 시각화
ggplot(df_heatmap[1:200, ], 
       aes(x = X, y = Y, fill = Median)) +
  geom_tile(color = "white") +  # 타일 형태로 색상 적용
  scale_fill_viridis_c(option = "plasma") +  # 색상 스케일 지정
  facet_wrap(~ Year, ncol = 5) +  # 연도별로 분할
  labs(title = "5×5 Grid - 90th Percentile of Precipitation by Year",
       x = "Grid X",
       y = "Grid Y",
       fill = "90th Percentile") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"))










# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run1[1,1,]  

# 상관계수 저장 벡터
cor_results <- numeric(365)  

# i < 365에 대해 연도별 데이터 추출 및 상관계수 계산
for (i in 0:364) {
  # 164년치 동일 날짜 데이터 추출 (마지막 연도 제외)
  series1 <- data_series[i + (0:163) * 365]      
  series2 <- data_series[i + 1 + (0:163) * 365*1]  
  
  # 길이가 같은지 확인
  if (length(series1) == length(series2)) {
    cor_results[i + 1] <- cor(series1, series2, use="complete.obs")
  } else {
    cor_results[i + 1] <- NA  # 길이가 다르면 NA 처리
  }
}

# 결과 확인
plot(0:364, cor_results, type="l", main="1년 차이 상관계수", 
     xlab="날짜 (i)", ylab="상관계수", col="blue")




# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run1[1,1,]

# 연도 개수 계산
num_years <- length(data_series) / 365  # 총 연도 수 (약 165년)

# 월별 데이터 저장 리스트
monthly_data <- vector("list", 12)

# 월별 데이터 추출 (각 연도의 동일한 월 데이터)
for (m in 1:12) {
  monthly_data[[m]] <- sapply(0:(num_years-1), function(y) {
    start_idx <- y * 365 + sum(c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)[m]) + 1
    end_idx <- start_idx + c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[m] - 1
    return(mean(data_series[start_idx:end_idx], na.rm = TRUE))
  })
}

# 월별 주기성 분석 (ACF 활용)
par(mfrow=c(4,3))  # 12개 그래프를 위한 배치 설정
for (m in 1:12) {
  acf(monthly_data[[m]], main=paste(m, "월 ACF 분석"), lag.max=50, ylim=c(-1,1))
}



















# 예제 데이터 (run1이 주어진 상태라고 가정)
# run1은 (5,5,60225) 형태의 다차원 배열

# 1. 5x5 지역의 합 벡터 계산
total_sum <- apply(run1, 3, sum)  # (60225,) 형태

# 2. (1,1) 위치의 강수량 데이터
target_series <- run1[1,1,]  # (60225,) 형태

# 3. 상관계수 계산
correlation <- cor(total_sum, target_series, method = "pearson")

# 결과 출력
print(correlation)

plot(total_sum, target_series)

plot(run1[5,5,], target_series)









# (1,1) 위치의 강수량 데이터
X <- run1[1,1,]  # (60225,) 형태

# 1차 차분 벡터 생성
diff1 <- diff(X, differences = 1)  # (60224,) 형태

# 한 시점 이동된 1차 차분 벡터 생성
diff2 <- diff1[-1]  # (60223,) 형태
diff1 <- diff1[-length(diff1)]  # 같은 길이로 맞춤 (60223,)

# 상관계수 계산
correlation <- cor(diff1, diff2, method = "pearson")

# 결과 출력
print(correlation)









hist(run1[1,1,])


# (1,1) 위치의 강수량 데이터
X_t <- run1[2,4,]  # (60225,)

# 시계열 길이
T <- length(X_t)

# X_{t-1}과 X_t 데이터 생성
X_t_minus1 <- X_t[1:(T-1)]  # 이전 시점 데이터
X_t_current <- X_t[2:T]      # 현재 시점 데이터

# 상관계수 계산
correlation <- cor(X_t_minus1, X_t_current, method = "pearson")

# 결과 출력
cat("X_t와 X_{t-1}의 상관계수:", correlation, "\n")



# 데이터 프레임 생성
df <- data.frame(X_t_minus1 = X_t_minus1, X_t_current = X_t_current)

# 밀도를 색으로 표현한 산점도 (2D binning)
ggplot(df, aes(x = X_t_minus1, y = X_t_current)) +
  geom_bin2d(bins = 50) +  # bins 개수를 조절하여 해상도 조정
  scale_fill_gradient(low = "blue", high = "red") +  # 농도를 파랑~빨강으로 표현
  labs(x = "X_{t-1}", y = "X_t", title = "X_{t-1} vs X_t (Density Plot)") +
  theme_minimal()




# (1,1) 위치의 강수량 데이터
X_t <- run1[1,1,]  # (60225,)

# 1. 시계열 그래프 (최근 1000일 데이터 확인)
ts_data <- ts(X_t)  # 시계열 데이터 변환
plot(ts_data[1:1000], type="l", col="blue", main="Time Series Plot (최근 1000일)",
     xlab="Time", ylab="Precipitation")

# 2. 자기상관함수(ACF) 분석
acf(X_t, lag.max=365, main="Autocorrelation Function (ACF)")  # 최대 1년 주기 확인

# 3. 주파수 분석 (FFT)
spectrum <- abs(fft(X_t))  # FFT 변환 후 절대값 계산
freq <- (0:(length(spectrum)-1)) / length(spectrum)  # 주파수 값 설정
plot(freq[1:1000], spectrum[1:1000], type="l", col="red",
     main="Frequency Spectrum (FFT)",
     xlab="Frequency", ylab="Magnitude")




# (1,1) 위치의 강수량 데이터
X_t <- run1[1,1,]  

# Fast Fourier Transform (FFT) 수행
spectrum <- abs(fft(X_t))  # FFT 변환 후 절대값 계산
N <- length(X_t)  # 데이터 길이

# 주파수 계산
freq <- (0:(N-1)) / N  # 정규화된 주파수
periods <- 1 / freq  # 주기에 해당하는 값 (0 제외)

# 데이터 프레임 생성 (주파수와 스펙트럼 값 저장)
fft_data <- data.frame(Frequency = freq, Period = periods, Magnitude = spectrum)

# 0 주기(무한대)를 제외하고 크기 순으로 정렬
fft_data <- fft_data[-1, ]  # 첫 번째 값(DC component) 제거
fft_data <- fft_data[order(-fft_data$Magnitude), ]  # 크기 순 정렬

# 상위 10개 주기 출력
library(dplyr)
top_periods <- fft_data %>%
  filter(Period < N/2) %>%  # 비현실적으로 큰 주기 제거
  head(10)  # 상위 10개 주기 선택

# 결과 출력
print(top_periods)












# (1,1) 위치의 강수량 데이터
X_t <- run1[1,1,]  
T <- length(X_t)

# 1. 데이터 프레임 생성 (최근 5000일만 사용)
df <- data.frame(Time = 1:5000, Precipitation = X_t[1:5000])

# 1️⃣ **원래 시계열 플롯**
ggplot(df, aes(x = Time, y = Precipitation)) +
  geom_line(color = "blue") +
  labs(title = "Original Time Series (First 5000 Days)",
       x = "Time", y = "Precipitation") +
  theme_minimal()

# 2️⃣ **700일 이동평균 플롯**
rolling_avg <- rollmean(X_t, k=365, fill=NA)  # 700일 이동평균 계산
df$Rolling_Avg <- rolling_avg[1:5000]  # 데이터프레임에 추가

ggplot(df[1:5000,], aes(x = Time)) +
  # geom_line(aes(y = Precipitation), color = "blue", alpha = 0.4) +  # 원래 데이터
  geom_line(aes(y = Rolling_Avg), color = "red", size = 1) +  # 이동평균
  labs(title = "700-day Moving Average",
       x = "Time", y = "Smoothed Precipitation") +
  theme_minimal()

# 3️⃣ **700일 주기 필터링 그래프 (FFT)**
spectrum <- fft(X_t)  # FFT 변환
freq <- (0:(T-1)) / T  # 주파수 계산
target_freq <- which.min(abs(1/freq - 700))  # 700일 주기에 해당하는 주파수 찾기

# 700일 주기만 남기고 필터링
filtered_spectrum <- rep(0+0i, length(spectrum))
filtered_spectrum[target_freq] <- spectrum[target_freq]  # 700일 주기만 남김
filtered_signal <- Re(fft(filtered_spectrum, inverse=TRUE) / length(spectrum))  # 역변환

# 필터링된 데이터프레임 생성
df$Filtered_700 <- filtered_signal[1:5000]

ggplot(df, aes(x = Time, y = Filtered_700)) +
  geom_line(color = "purple", size = 1) +
  labs(title = "Filtered 700-day Cycle",
       x = "Time", y = "Precipitation (Filtered)") +
  theme_minimal()





# (1,1) 위치의 강수량 데이터
X_t <- run1[3,1,]  # (60225,)

# 시계열 길이
T <- length(X_t)

# X_{t-1}과 X_t 데이터 생성
X_t_minus1 <- X_t[1:(T-1)]  # 이전 시점 데이터
X_t_current <- X_t[2:T]      # 현재 시점 데이터

# 데이터 프레임 생성
df <- data.frame(X_t_minus1 = X_t_minus1, X_t_current = X_t_current)

# 밀도를 색으로 표현한 산점도 (2D binning)
ggplot(df, aes(x = X_t_minus1, y = X_t_current)) +
  geom_bin2d(bins = 50) +  # bins 개수를 조절하여 해상도 조정
  scale_fill_gradient(low = "blue", high = "red") +  # 농도를 파랑~빨강으로 표현
  labs(x = "X_{t-1}", y = "X_t", title = "X_{t-1} vs X_t (Density Plot)") +
  theme_minimal()




# 빈 리스트 생성 (각 위치별 ggplot 저장)
plot_list <- list()
# 5x5 격자에서 각각의 위치에 대해 그래프 생성
for (i in 1:5) {
  for (j in 1:5) {
    # (i,j) 위치의 강수량 데이터
    X_t <- run4[i,j,]  
    T <- length(X_t)
    
    # X_{t-1}과 X_t 데이터 생성
    X_t_minus1 <- X_t[1:(T-1)]
    X_t_current <- X_t[2:T]
    
    # 데이터 프레임 생성
    df <- data.frame(X_t_minus1 = X_t_minus1, X_t_current = X_t_current)
    
    # 밀도를 색으로 표현한 산점도 (범례 및 타이틀 제거)
    p <- ggplot(df, aes(x = X_t_minus1, y = X_t_current)) +
      geom_bin2d(bins = 50) +  
      scale_fill_gradient(low = "blue", high = "red", guide = "none") +  # 범례 제거
      labs(x = NULL, y = NULL) +  # 축 라벨 제거
      theme_minimal() +
      theme(legend.position = "none",  # 범례 제거
            plot.title = element_blank(),  # 타이틀 제거
            axis.text = element_blank(),  # 축 눈금 제거
            axis.ticks = element_blank())  # 축 눈금 제거
    
    # 리스트에 저장
    plot_list[[length(plot_list) + 1]] <- p
  }
}

x_label <- textGrob(expression(X[t-1]), gp = gpar(fontsize = 14, fontface = "bold"))
y_label <- textGrob(expression(X[t]), rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
title_label = textGrob("run4", gp = gpar(fontsize = 14, fontface = "bold"))

# 5x5 격자로 플롯 정렬
grid.arrange(
  arrangeGrob(grobs = plot_list, ncol = 5, nrow = 5),
  bottom = x_label,  # X축 제목
  left = y_label,     # Y축 제목
  top = title_label
)


X_t_minus1[which( X_t_current >3)]
X_t_current[which( X_t_current >3)]

X_t_minus1[which( X_t_minus1 >3)]
X_t_current[which( X_t_minus1 >3)] %>% round(3)




















# 특정 지역의 시계열 데이터를 저장할 행렬 (5×5 지역)
cor_matrix <- matrix(0, 25, 25)  # 25개 지역 간 상관행렬

# 지역 간 시계열 상관관계 계산
for (i in 1:5) {
  for (j in 1:5) {
    for (k in 1:5) {
      for (l in 1:5) {
        # 1D 인덱스로 변환
        idx1 <- (i - 1) * 5 + j
        idx2 <- (k - 1) * 5 + l
        
        # 피어슨 상관계수 계산
        cor_matrix[idx1, idx2] <- cor(run1[i,j,], run1[k,l,], use="complete.obs")
      }
    }
  }
}

# 결과 시각화 (Heatmap)
library(ggplot2)
library(reshape2)

# 데이터 프레임 변환
cor_df <- melt(cor_matrix)
colnames(cor_df) <- c("Region1", "Region2", "Correlation")

# 히트맵 그리기
ggplot(cor_df, aes(x=Region1, y=Region2, fill=Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal() +
  labs(title="지역 간 시계열 상관행렬", x="지역1", y="지역2")












# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run1[1,1,]

# X_{t-1}과 변화량 (X_t - X_{t-1}) 계산
X_t1 <- data_series[1:(length(data_series)-1)]
X_diff <- diff(data_series)  # X_t - X_{t-1}

# 데이터프레임 생성
df <- data.frame(X_t1, X_diff)

# X_{t-1} 값을 구간(bin)으로 나눠 평균 변화량 계산
library(dplyr)
df_summary <- df %>%
  mutate(X_bin = cut(X_t1, breaks=20)) %>%  # 20개 구간으로 나눔
  group_by(X_bin) %>%
  summarise(mean_diff = mean(X_diff, na.rm=TRUE))

# 바 플롯 그리기
library(ggplot2)
ggplot(df_summary, aes(x=X_bin, y=mean_diff)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="X_{t-1}에 따른 변화량 (X_t - X_{t-1})", 
       x="X_{t-1} (구간)", y="X_t - X_{t-1} 평균 변화량") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





















library(dplyr)

# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run1[1,1,]

# X_{t-1}과 X_t 생성
X_t1 <- data_series[1:(length(data_series)-1)]
X_t <- data_series[2:length(data_series)]

# 데이터프레임 생성
df <- data.frame(X_t1, X_t)

# X_{t-1}을 분위수(quantile) 기준으로 4개 그룹으로 나눔
df <- df %>%
  mutate(X_bin = ntile(X_t1, 4))  # 4분위(Quartile)로 구분

# 각 구간별로 AR(1) 모델 적합
ar_results <- df %>%
  group_by(X_bin) %>%
  summarise(beta1 = coef(lm(X_t ~ X_t1, data = cur_data()))[2],
            p_value = summary(lm(X_t ~ X_t1, data = cur_data()))$coefficients[2,4])

# 결과 출력
print(ar_results)

# 계수 시각화
library(ggplot2)
ggplot(ar_results, aes(x=factor(X_bin), y=beta1)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="X_{t-1} 크기에 따른 AR(1) 계수 변화", x="X_{t-1} 구간", y="AR(1) 계수")

