library(dplyr)
library(zoo)
library(ggplot2)
library(reshape2)
library(lubridate)
library(tidyr)
library(moments)  # 첨도 계산을 위한 패키지

load("data//run1.RData")
load("data//run2.RData")
load("data//run3.RData")
load("data//run4.RData")


# 1. 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)
yearSet = year(dates) %>% unique()
years = year(dates)

# 2. run1[1,1,]의 데이터를 y값으로 사용
y_values <- run1[1,1,]

# 3. 90일 이동 중앙값(Moving Median) 계산
y_moving_median <- rollapply(y_values, width = 90, FUN = median, fill = NA, align = "center")

# 4. 90일 이동 90퍼센타일(Moving 90th Percentile) 계산
y_moving_90th <- rollapply(y_values, width = 90, FUN = function(x) quantile(x, 0.9), fill = NA, align = "center")

# 5. 데이터 프레임 생성
df <- data.frame(date = dates, median_value = y_moving_median, percentile_90 = y_moving_90th , value = y_values)

# 6. 9월 1일마다 vline 추가 (매년 9월 1일을 찾기)
september_1_dates <- seq(from = as.Date("1850-09-01"), to = as.Date("2014-09-01"), by = "year")

# 7. 시각화
ggplot(df[1:4000,], aes(x = date)) +
  # geom_line(aes(y=value))+
  geom_line(aes(y = median_value), color = "blue", size = 1, alpha = 0.8) +  # 90일 이동 중앙값
  geom_line(aes(y = percentile_90), color = "green", size = 1, alpha = 0.8) +  # 90일 이동 90퍼센타일
  # geom_vline(xintercept = as.numeric(september_1_dates), color = "red", linetype = "dashed") +  # 9월 1일마다 빨간 점선
  labs(title = "90-Day Moving Median & 90th Percentile Time Series",
       x = "Date", y = "Value",
       color = "Legend") +
  scale_color_manual(values = c("blue" = "Median", "green" = "90th Percentile")) +
  theme_minimal()




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










