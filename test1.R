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
library(plotly)
library(gstat)
library(sp)
library(spacetime)
library(sf)
library(geosphere)
library(evd)
library(evir)   # gpd() 함수 제공
library(fields)
library(lmtest)

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
  DayOfYear = yday(dates)  # 연도의 몇 번째 날인지
)

# 변환된 시간 계산 (연도를 소수점으로 변환)
df_time <- df_time %>%
  mutate(Time_Normalized = (Year-1850) + (DayOfYear / 60225))




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







# 5x5 셀의 평균 강수량 계산
# 5x5 셀의 평균 강수량 계산
mean_precipitation <- apply(run1, c(1, 2), mean)

# X, Y 좌표 생성
x_vals <- seq(1, 5, 1)
y_vals <- seq(1, 5, 1)

# plotly의 surface 플롯을 위한 올바른 데이터 형태로 변환
fig <- plot_ly(
  x = x_vals,
  y = y_vals,
  z = ~mean_precipitation,
  type = "surface"
)

# 그래프 레이아웃 설정
fig <- fig %>% layout(
  title = "3D Mean Precipitation Distribution",
  scene = list(
    xaxis = list(title = "X-axis (Grid)"),
    yaxis = list(title = "Y-axis (Grid)"),
    zaxis = list(title = "Mean Precipitation")
  )
)

# 그래프 출력
fig


#####################################################
# 5x5 평균 강수량 계산 함수
calculate_mean_precip <- function(run_data, label) {
  mean_precip <- apply(run_data, c(1, 2), mean)
  df <- melt(mean_precip)
  colnames(df) <- c("X", "Y", "MeanPrecip")
  df$Run <- label
  return(df)
}

# 각 Run 데이터 변환
df1 <- calculate_mean_precip(run1, "Run1")
df2 <- calculate_mean_precip(run2, "Run2")
df3 <- calculate_mean_precip(run3, "Run3")
df4 <- calculate_mean_precip(run4, "Run4")

# 전체 데이터에서 최소, 최대값 구하기
all_data <- rbind(df1, df2, df3, df4)
min_val <- min(all_data$MeanPrecip)
max_val <- max(all_data$MeanPrecip)

# 공통 스타일을 적용한 히트맵 생성 함수 (동일한 컬러 스케일 범위 설정)
create_heatmap <- function(df, title) {
  ggplot(df, aes(x = X, y = Y, fill = MeanPrecip)) +
    geom_tile() +
    scale_fill_gradient(low = "grey", high = "red", limits = c(min_val, max_val)) +  # 동일한 범위 설정
    labs(title = title, x = "X-axis (Grid)", y = "Y-axis (Grid)", fill = "Mean Precip") +
    theme_minimal()
}

# 4개 히트맵 생성 (컬러 스케일 동일)
plot1 <- create_heatmap(df1, "Run1 - Mean Precipitation")
plot2 <- create_heatmap(df2, "Run2 - Mean Precipitation")
plot3 <- create_heatmap(df3, "Run3 - Mean Precipitation")
plot4 <- create_heatmap(df4, "Run4 - Mean Precipitation")

# 4개 그래프를 한 화면에 배치
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)




#########################

# 하위 X% 설정
lower_percentile <- 1 # 예: 하위 10% 선택

# 하위 X% 평균 강수량 계산 함수
calculate_lower_percentile_mean <- function(run_data, label, percentile) {
  lower_mean_precip <- matrix(NA, nrow = 5, ncol = 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      # 특정 셀의 강수량 데이터 추출
      values <- run_data[i, j, ]
      
      # 하위 X% 값 선택
      threshold <- quantile(values, probs = percentile / 100, na.rm = TRUE)
      lower_values <- values[values <= threshold]
      
      # 하위 X% 평균 계산
      lower_mean_precip[i, j] <- mean(lower_values, na.rm = TRUE)
    }
  }
  
  # 데이터 프레임 변환
  df <- melt(lower_mean_precip)
  colnames(df) <- c("X", "Y", "MeanPrecip")
  df$Run <- label
  return(df)
}

# 각 Run 데이터 변환 (하위 X% 평균 강수량)
df1 <- calculate_lower_percentile_mean(run1, "Run1", lower_percentile)
df2 <- calculate_lower_percentile_mean(run2, "Run2", lower_percentile)
df3 <- calculate_lower_percentile_mean(run3, "Run3", lower_percentile)
df4 <- calculate_lower_percentile_mean(run4, "Run4", lower_percentile)

# 전체 데이터에서 최소, 최대값 구하기 (컬러 스케일 통일)
all_data <- rbind(df1, df2, df3, df4)
min_val <- min(all_data$MeanPrecip, na.rm = TRUE)
max_val <- max(all_data$MeanPrecip, na.rm = TRUE)

# 공통 스타일을 적용한 히트맵 생성 함수 (동일한 컬러 스케일 범위 설정)
create_heatmap <- function(df, title) {
  ggplot(df, aes(x = X, y = Y, fill = MeanPrecip)) +
    geom_tile() +
    scale_fill_gradient(low = "grey", high = "red", limits = c(min_val, max_val)) +  # 동일한 범위 설정
    labs(title = title, x = "X-axis (Grid)", y = "Y-axis (Grid)", fill = "Mean Precip") +
    theme_minimal()
}

# 4개 히트맵 생성 (컬러 스케일 동일)
plot1 <- create_heatmap(df1, paste0("Run1 - Lower ", lower_percentile, "% Mean Precip"))
plot2 <- create_heatmap(df2, paste0("Run2 - Lower ", lower_percentile, "% Mean Precip"))
plot3 <- create_heatmap(df3, paste0("Run3 - Lower ", lower_percentile, "% Mean Precip"))
plot4 <- create_heatmap(df4, paste0("Run4 - Lower ", lower_percentile, "% Mean Precip"))

# 4개 그래프를 한 화면에 배치
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)





# 분석할 특정 위치 선택 (예: (1,1))
i <- 1
j <- 1

# run1[i, j, ]에서 1차원 벡터 추출
rainfall_data <- run1[i, j, ]

# 데이터 정보 설정
days_per_year <- 365
years <- length(rainfall_data) / days_per_year  # 전체 연수 (165년)

# 날짜 벡터 생성
day_of_year <- rep(1:days_per_year, years)
year <- rep(1850:(1850 + years - 1), each=days_per_year)

# 데이터 프레임 생성
df <- data.frame(Day = factor(day_of_year), Year = year, Rainfall = rainfall_data)

# 박스플롯 생성
ggplot(df, aes(x=Day, y=Rainfall)) +
  geom_boxplot(outlier.shape=NA) +  # 이상치 제거
  labs(title="Daily Precipitation Distribution (1850-2014)", 
       x="Day of the Year", y="Precipitation (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) 

####################################################################################





# 데이터 정보 설정
days_per_year <- 365
years <- 165  # 1850년~2014년 (총 165년)
grid_size <- 5  # 5x5 지역

# 분위수 데이터를 저장할 리스트
percentile_list <- list()

# 5x5 지역에 대해 반복문 수행
for (i in 1:grid_size) {
  for (j in 1:grid_size) {
    # 각 위치의 강수량 데이터 가져오기
    rainfall_data <- run1[i, j, ]
    
    # 날짜 벡터 생성
    day_of_year <- rep(1:days_per_year, years)
    
    # 데이터 프레임 생성
    df <- data.frame(Cell = paste0("(", i, ",", j, ")"), 
                     Day = day_of_year, 
                     Rainfall = rainfall_data)
    
    # 분위수(10%, 25%, 50%, 75%, 90%) 계산
    percentiles <- df %>%
      group_by(Day) %>%
      summarise(p10 = quantile(Rainfall, 0.10, na.rm=TRUE),
                p25 = quantile(Rainfall, 0.25, na.rm=TRUE),
                p50 = quantile(Rainfall, 0.50, na.rm=TRUE),
                p75 = quantile(Rainfall, 0.75, na.rm=TRUE),
                p90 = quantile(Rainfall, 0.90, na.rm=TRUE)) %>%
      mutate(Cell = paste0("(", i, ",", j, ")"))  # 셀 위치 추가
    
    # 리스트에 저장
    percentile_list[[paste0(i, "_", j)]] <- percentiles
  }
}

# 모든 데이터를 하나의 데이터 프레임으로 결합
percentile_df <- bind_rows(percentile_list)

# X축 눈금 (월별 대표 일자)
month_ticks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# 5x5 지역에 대해 분위수 그래프 생성 (facet_wrap 사용)
ggplot(percentile_df, aes(x=Day)) +
  geom_line(aes(y=p10, color="10%"), linewidth=0.6) +
  geom_line(aes(y=p25, color="25%"), linewidth=0.6) +
  geom_line(aes(y=p50, color="50% (Median)"), linewidth=0.8) +
  geom_line(aes(y=p75, color="75%"), linewidth=0.6) +
  geom_line(aes(y=p90, color="90%"), linewidth=0.6) +
  labs(title="Precipitation Percentiles by Day (1850-2014) for 5x5 Grid",
       x="Day of the Year", y="Precipitation (mm)", color="Percentile") +
  theme_minimal() +
  scale_x_continuous(breaks=month_ticks, labels=month_labels) +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~Cell, ncol=5)










# 5x5 지역의 가능한 좌표 생성
coords <- expand.grid(X = 1:5, Y = 1:5)

# 날짜 벡터 생성 (1850-2014년)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)

# 15일 단위 블록 생성
time_blocks <- rep(1:(length(dates) / 15), each=15, length.out=length(dates))

# 시간 블록을 공간 데이터와 결합
rainfall_data <- expand.grid(X = 1:5, Y = 1:5, Date = dates) %>%
  mutate(Block = rep(time_blocks, times = 5 * 5))  # 각 위치(X, Y)에 대해 반복 적용

# 강수량 데이터를 올바르게 할당
rainfall_data$Rainfall <- apply(rainfall_data, 1, function(row) {
  X_idx <- as.integer(row["X"])
  Y_idx <- as.integer(row["Y"])
  date_idx <- which(dates == as.Date(row["Date"]))  # 날짜 인덱스 찾기
  if (length(date_idx) == 1 && X_idx >= 1 && X_idx <= 5 && Y_idx >= 1 && Y_idx <= 5) {
    return(run1[X_idx, Y_idx, date_idx])
  } else {
    return(NA)  # 범위를 벗어나면 NA
  }
})

# 가능한 모든 위치 쌍 계산
location_pairs <- expand.grid(1:nrow(coords), 1:nrow(coords))
colnames(location_pairs) <- c("point1", "point2")

# 유클리드 거리 계산 (공간적 거리)
location_pairs <- location_pairs %>%
  rowwise() %>%
  mutate(
    X1 = coords$X[point1], Y1 = coords$Y[point1],
    X2 = coords$X[point2], Y2 = coords$Y[point2],
    distance = sqrt((X2 - X1)^2 + (Y2 - Y1)^2)
  ) %>%
  ungroup()

# 강수량 데이터와 위치 데이터를 조인하여 거리별 분석
variogram_data <- rainfall_data %>%
  inner_join(rainfall_data, by = c("Block"), suffix = c("_1", "_2")) %>%
  inner_join(location_pairs, by = c("X_1" = "X1", "Y_1" = "Y1", "X_2" = "X2", "Y_2" = "Y2"))

# 반분산(Semivariance) 계산
variogram_data <- variogram_data %>%
  mutate(semivariance = (Rainfall_1 - Rainfall_2)^2 / 2) %>%
  group_by(Block, distance) %>%
  summarise(semivariance = mean(semivariance, na.rm=TRUE), .groups="drop")

# 15일 단위 Variogram 평균 계산
final_variogram <- variogram_data %>%
  group_by(distance) %>%
  summarise(semivariance = mean(semivariance, na.rm=TRUE), .groups="drop")

# ggplot을 이용한 Variogram 시각화
ggplot(final_variogram, aes(x=distance, y=semivariance)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="loess", se=FALSE, color="red") +
  labs(title="Averaged Empirical Spatial Variogram (15-day blocks)",
       x="Spatial Lag Distance",
       y="Semivariance") +
  theme_minimal()











###################################################

library(evir)   # gpd() 함수 제공
data = as.vector(run1) 


# 여러 개의 threshold 값 설정 (상위 90% ~ 99% 분위수)
thresholds <- quantile(data, seq(0.90, 0.99, by = 0.01))

# threshold 별로 GPD 적합 후 모수 추정값 저장
shape_values <- c()
scale_values <- c()

for (thresh in thresholds) {
  fit <- gpd(data, thresh)  # GPD 적합
  shape_values <- c(shape_values, fit$par.ests["xi"])  # Shape (ξ)
  scale_values <- c(scale_values, fit$par.ests["beta"]) # Scale (σ)
}

# 결과를 데이터 프레임으로 정리
results <- data.frame(Threshold = thresholds, Shape = shape_values, Scale = scale_values)

# Shape 변화 시각화
p1 <- ggplot(results, aes(x = Threshold, y = Shape)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue") +
  labs(title = "Threshold 변화에 따른 Shape (ξ) 추정값",
       x = "Threshold",
       y = "Shape (ξ)") +
  theme_minimal()

# Scale 변화 시각화
p2 <- ggplot(results, aes(x = Threshold, y = Scale)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "red") +
  labs(title = "Threshold 변화에 따른 Scale (σ) 추정값",
       x = "Threshold",
       y = "Scale (σ)") +
  theme_minimal()

# 두 그래프 나란히 배치
grid.arrange(p1, p2, ncol = 2)


#######################################

# 데이터 준비
data <- as.vector(run1)  # (5,5,60225) 데이터를 1차원 벡터로 변환
total_days <- length(data)  # 전체 일수 (60225)
years <- total_days / 365  # 연도 개수

# 계절 정의 (봄: 3~5월, 여름: 6~8월, 가을: 9~11월, 겨울: 12~2월)
season_names <- c("Spring", "Summer", "Autumn", "Winter")
season_indices <- list(
  Spring = rep(c(rep(TRUE, 59), rep(FALSE, 274)), years),
  Summer = rep(c(rep(FALSE, 59), rep(TRUE, 92), rep(FALSE, 214)), years),
  Autumn = rep(c(rep(FALSE, 151), rep(TRUE, 91), rep(FALSE, 123)), years),
  Winter = rep(c(rep(FALSE, 242), rep(TRUE, 31), rep(FALSE, 59), rep(TRUE, 31)), years)
)

# 분석할 threshold 비율 (상위 90% ~ 99% 분위수)
quantiles <- seq(0.90, 0.99, by = 0.01)

# 결과 저장 리스트
seasonal_results <- list()

# 계절별 분석 수행
for (season in season_names) {
  seasonal_data <- data[season_indices[[season]]]  # 계절별 데이터 추출
  
  shape_values <- c()
  scale_values <- c()
  thresholds <- quantile(seasonal_data, quantiles)  # 계절별 threshold 계산
  
  # 각 threshold에 대해 GPD 적합
  for (thresh in thresholds) {
    fit <- gpd(seasonal_data, thresh)  # GPD 적합
    shape_values <- c(shape_values, fit$par.ests["xi"])   # Shape (ξ)
    scale_values <- c(scale_values, fit$par.ests["beta"]) # Scale (σ)
  }
  
  # 결과 저장
  seasonal_results[[season]] <- data.frame(
    Threshold = thresholds,
    Shape = shape_values,
    Scale = scale_values,
    Season = season
  )
}

# 전체 데이터프레임 결합
final_results <- do.call(rbind, seasonal_results)

# 개별 그래프 생성
p1 <- ggplot(final_results, aes(x = Threshold, y = Shape, color = Season)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Threshold 변화에 따른 Shape (ξ) 추정값",
       x = "Threshold",
       y = "Shape (ξ)",
       color = "Season") +
  theme_minimal()

p2 <- ggplot(final_results, aes(x = Threshold, y = Scale, color = Season)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Threshold 변화에 따른 Scale (σ) 추정값",
       x = "Threshold",
       y = "Scale (σ)",
       color = "Season") +
  theme_minimal()

# 그래프 나란히 출력
grid.arrange(p1, p2, ncol = 2)




#########################################


# 데이터 준비
data <- as.vector(run1)  # (5,5,60225) 데이터를 1차원 벡터로 변환
total_days <- length(data)  # 전체 일수 (60225)
years <- total_days / 365  # 연도 개수

# 계절 정의 (봄: 3~5월, 여름: 6~8월, 가을: 9~11월, 겨울: 12~2월)
season_names <- c("Spring", "Summer", "Autumn", "Winter")
season_indices <- list(
  Spring = rep(c(rep(TRUE, 59), rep(FALSE, 274)), years),
  Summer = rep(c(rep(FALSE, 59), rep(TRUE, 92), rep(FALSE, 214)), years),
  Autumn = rep(c(rep(FALSE, 151), rep(TRUE, 91), rep(FALSE, 123)), years),
  Winter = rep(c(rep(FALSE, 242), rep(TRUE, 31), rep(FALSE, 59), rep(TRUE, 31)), years)
)

# 계절별 강수량 데이터 추출
season_data <- data.frame(Precipitation = numeric(0), Season = character(0))

for (season in season_names) {
  season_values <- data[season_indices[[season]]]  # 해당 계절의 강수량 데이터 추출
  temp_df <- data.frame(Precipitation = season_values, Season = season)
  season_data <- rbind(season_data, temp_df)  # 데이터 병합
}

# 계절별 분산 계산
season_variance <- season_data %>%
  group_by(Season) %>%
  summarise(Variance = var(Precipitation, na.rm = TRUE))

# 결과 출력
print(season_variance)

# 박스플롯으로 계절별 강수량 분산 비교
ggplot(season_data, aes(x = Season, y = Precipitation, fill = Season)) +
  geom_boxplot() +
  labs(title = "계절별 강수량 분포", x = "계절", y = "강수량") +
  theme_minimal()

# Bartlett's Test (정규성 가정)
bartlett_test <- bartlett.test(Precipitation ~ Season, data = season_data)
print(bartlett_test)


#######################################################################








# 365일 주기 기반 푸리에 회귀 모델 생성
fourier_terms <- function(t, K = 3, period = 365) {
  terms <- data.frame(t = t)
  for (k in 1:K) {
    terms[[paste0("sin", k)]] <- sin(2 * pi * k * t / period)
    terms[[paste0("cos", k)]] <- cos(2 * pi * k * t / period)
  }
  return(terms)
}

# 3차 푸리에 근사 (K = 3)
K <- 2
fourier_data <- fourier_terms(t, K)

# 선형 회귀 모델 적합
model <- lm(X_t ~ ., data = fourier_data)

# 계절성 성분 S(t) 예측
S_t <- predict(model, newdata = fourier_data)

# 결과 확인
print(length(S_t))  # 36500일의 계절성 성분






####################################

df <- df %>%
  mutate(
    season = case_when(
      month(date) %in% c(3,4,5,6,7,8) ~ "Warm Season",  # 봄+여름
      month(date) %in% c(9,10,11,12,1,2) ~ "Cold Season" # 가을+겨울
    )
  )


anova_result <- aov(precip ~ season, data = df)
summary(anova_result)

kruskal_result <- kruskal.test(precip ~ season, data = df)
kruskal_result



###################################################


# 데이터를 long format으로 변환
df <- expand.grid(i = 1:5, j = 1:5, date = dates) %>%
  mutate(
    doy = yday(date),  # 연중일 (1~365)
    precip = as.vector(run1)  # 강수량 데이터
  )

# 연중일(DOY)별 분위수 계산 (5x5 셀 전체에 대해 집계)
quantile_df <- df %>%
  group_by(doy) %>%
  summarise(
    q10 = quantile(precip, 0.10),
    q50 = quantile(precip, 0.50),
    mean = mean(precip),
    q90 = quantile(precip, 0.90),
    q99 = quantile(precip, 0.99)
  ) %>%
  pivot_longer(cols = c(q10, q50, mean, q90, q99), names_to = "quantile", values_to = "value")

# 시각화
ggplot(quantile_df, aes(x = doy, y = value, color = quantile)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(0, 365, 30)) +  # X축을 월별로 보기 쉽게 조정
  labs(
    title = "Annual Precipitation Quantiles (5×5 Grid)",
    x = "Day of Year (DOY)",
    y = "Precipitation",
    color = "Quantile"
  ) +
  theme_minimal()










