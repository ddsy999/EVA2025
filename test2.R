# 5x5 블록 단위로 합산하는 함수
block_sum <- function(arr, block_size = 5) {
  dim_x <- dim(arr)[1]
  dim_y <- dim(arr)[2]
  dim_t <- dim(arr)[3]
  
  # 블록 개수 계산
  new_x <- dim_x %/% block_size
  new_y <- dim_y %/% block_size
  
  # 블록 합 벡터 생성
  result <- numeric(new_x * new_y * dim_t)
  index <- 1
  
  for (t in 1:dim_t) {
    for (i in seq(1, dim_x, by = block_size)) {
      for (j in seq(1, dim_y, by = block_size)) {
        # 5x5 블록의 합을 계산
        block_values <- arr[i:(i+block_size-1), j:(j+block_size-1), t]
        result[index] <- sum(block_values, na.rm = TRUE)
        index <- index + 1
      }
    }
  }
  
  return(result)
}

# 벡터 변환 (길이: 60225)
run1_vector <- block_sum(run1, block_size = 5)
length(run1_vector)  # 60225 확인

##############################################################


# 데이터 길이 설정
N <- length(run1_vector)  # 총 데이터 개수 (60225)
days_per_year <- 365
years_per_group <- 10
days_per_group <- years_per_group * days_per_year  # 10년치 (3650일)

# 사용할 데이터 길이 (3650일씩 나누어 떨어지는 가장 가까운 값)
N_truncated <- floor(N / days_per_group) * days_per_group  # 58400개 (16개 그룹)

# 데이터를 3650일 단위로 reshape
run1_trimmed <- run1_vector[1:N_truncated]  # 마지막 625개 데이터 삭제
run1_matrix <- matrix(run1_trimmed, ncol = days_per_group, byrow = TRUE)

# 10년 단위 평균 계산 (세로 방향 평균)
ten_year_avg <- apply(run1_matrix, 2, mean, na.rm = TRUE)  # 3650개 데이터 생성

# 결과 확인
print(length(ten_year_avg))  # 3650개여야 함





# FFT 변환 수행
fft_result <- fft(ten_year_avg)
power_spectrum <- Mod(fft_result)^2  # 진폭의 제곱 (강도)
n <- length(ten_year_avg)  # 데이터 길이 (3650)

# 주파수 및 주기 계산
freqs <- (0:(n-1)) / n  # 주파수 값
periods <- 1 / freqs    # 주기 (일 단위)
valid_idx <- 2:(n/2)    # DC 성분(0번 주파수) 제거, Nyquist 대칭 부분 제외

# 데이터 프레임 생성
fft_df <- data.frame(Period = periods[valid_idx], Power = power_spectrum[valid_idx])

# 로그 스케일을 적용하여 FFT 분석 결과 시각화
ggplot(fft_df, aes(x = Period, y = Power)) +
  geom_point(color = "blue", size = 2) +                     # 점 그래프
  geom_line(linetype = "dashed", color = "blue", alpha = 0.7) +  # 선 연결
  scale_x_log10() +  # 로그 스케일 적용 (주기 간 차이 큼)
  labs(title = "FFT 분석 결과 - ten_year_avg 주기별 강도",
       x = "주기 (일)", 
       y = "강도") +
  theme_minimal()

# 강한 주기 3개 찾기
top_periods <- fft_df %>% arrange(desc(Power)) %>% slice(1:3)

# 상위 주기 출력
print(top_periods)

#################


# FFT 수행 (이미 계산한 fft_result 사용)
n <- length(ten_year_avg)  # 데이터 길이
fft_result <- fft(ten_year_avg)

# 강한 주기 3개
top_periods <- c(365.0, 182.5, 121.67)

# 해당 주기의 FFT 인덱스 찾기
top_indices <- round(n / top_periods)

# 진폭(A_k) 및 위상(phi_k) 계산
A_k <- 2 * Mod(fft_result[top_indices]) / n  # 진폭 정규화
phi_k <- Arg(fft_result[top_indices])  # 위상 (라디안)

# 절편값 (기본 평균값) 설정
beta_0 <- mean(ten_year_avg)

# 시간 벡터
t <- 1:n

# S_t 생성 (절편 포함한 사인 함수 모델)
S_t <- beta_0 + rowSums(sapply(1:length(top_periods), function(k) {
  A_k[k] * sin(2 * pi * t / top_periods[k] + phi_k[k])
}))

# 비교를 위해 데이터프레임 생성
df <- data.frame(
  t = 1:n,
  Actual = ten_year_avg,  # 실제 데이터
  Estimated_S_t = S_t  # 추정된 계절성 성분
)

# 시각화
ggplot(df, aes(x = t)) +
  geom_line(aes(y = Actual, color = "Actual Data")) +
  geom_line(aes(y = Estimated_S_t, color = "Estimated S(t) with Beta_0"), linetype = "dashed") +
  labs(title = "Comparison: Actual Data vs. Estimated S(t) with Intercept (Beta_0)",
       x = "Day", y = "Value") +
  scale_color_manual(values = c("Actual Data" = "blue", "Estimated S(t) with Beta_0" = "red")) +
  theme_minimal()



##################
# 데이터 길이 및 reshape
N <- length(run1_vector)  # 전체 데이터 길이 (36500일)
days_per_year <- 365
years <- N / days_per_year  # 총 100년

# 데이터 100년치를 365일 단위로 reshape
run1_matrix <- matrix(run1_vector, ncol = days_per_year, byrow = TRUE)

# 각 날짜별 평균 계산 (세로 방향 mean)
daily_mean <- colMeans(run1_matrix, na.rm = TRUE)  # 365개 생성

# S_t도 365일 단위로 reshape 후 평균 계산
S_t_reshaped <- matrix(S_t, ncol = days_per_year, byrow = TRUE)
S_t_avg <- colMeans(S_t_reshaped, na.rm = TRUE)

# 데이터프레임 생성
df <- data.frame(
  day = 1:days_per_year,
  Daily_Mean = daily_mean,  # 각 날짜의 평균
  S_t_Avg = S_t_avg  # S_t 평균값
)

# 비교 시각화
ggplot(df, aes(x = day)) +
  geom_line(aes(y = Daily_Mean, color = "Daily Mean (Actual)")) +
  geom_line(aes(y = S_t_Avg, color = "Estimated S(t)"), linetype = "dashed") +
  labs(title = "Comparison: Daily Mean vs. Estimated S(t)",
       x = "Day of Year", y = "Value") +
  scale_color_manual(values = c("Daily Mean (Actual)" = "blue", "Estimated S(t)" = "red")) +
  theme_minimal()




#######################


library(signal)
library(ggplot2)

# FFT 수행
n <- length(ten_year_avg)
fft_result <- fft(ten_year_avg)

# 상위 10개 강한 주기 찾기
valid_range <- 2:(n/2)
top_indices <- order(Mod(fft_result[valid_range])^2, decreasing = TRUE)[1:3] + 1  # 상위 10개 선택

# 해당 주기의 주파수 변환
top_periods <- n / top_indices

# 진폭(A_k) 및 위상(phi_k) 계산
A_k <- 2 * Mod(fft_result[top_indices]) / n  # 진폭 정규화
phi_k <- Arg(fft_result[top_indices])  # 위상 (라디안)

# 절편값 (기본 평균값)
beta_0 <- mean(ten_year_avg)

# 시간 벡터
t <- 1:n

# S_t 생성 (10개의 주기 사용)
S_t <- beta_0 + rowSums(sapply(1:length(top_periods), function(k) {
  A_k[k] * sin(2 * pi * t / top_periods[k] + phi_k[k])
}))

# 비교를 위해 데이터프레임 생성
df <- data.frame(
  t = 1:n,
  Actual = ten_year_avg,
  Estimated_S_t = S_t
)

# 시각화
ggplot(df[1:365*2,], aes(x = t)) +
  geom_line(aes(y = Actual, color = "Actual Data")) +
  geom_line(aes(y = Estimated_S_t, color = "Improved S(t)"), linetype = "dashed") +
  labs(title = "Comparison: Actual Data vs. Improved S(t)",
       x = "Day", y = "Value") +
  scale_color_manual(values = c("Actual Data" = "blue", "Improved S(t)" = "red")) +
  theme_minimal()

#############################

###########################################


df_ = data.frame(
  t = 1:365,
  Actual = daily_mean,
  Estimated_S_t = S_t[1:365]
)
# 시각화
ggplot(df_, aes(x = t)) +
  geom_line(aes(y = Actual, color = "Actual Data")) +
  geom_line(aes(y = Estimated_S_t, color = "Improved S(t)"), linetype = "dashed") +
  labs(title = "Comparison: Actual Data vs. Improved S(t)",
       x = "Day", y = "Value") +
  scale_color_manual(values = c("Actual Data" = "blue", "Improved S(t)" = "red")) +
  theme_minimal()






#############################
# 벡터 변환된 데이터 (5x5 블록 합산 후 60225개)
df <- data.frame(
  date = dates,
  doy = yday(dates),  # 연중일 (1~365)
  precip = run1_vector  # 5×5 블록을 합산한 강수량 값
)

# 연중일(DOY)별 분위수 계산 (전체 데이터에 대해 집계)
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
ggplot(quantile_df %>% filter(quantile =="mean"), aes(x = doy, y = value, color = quantile)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(0, 365, 30)) +  # X축을 월별로 보기 쉽게 조정
  labs(
    title = "Annual Precipitation Quantiles (After 5×5 Block Sum)",
    x = "Day of Year (DOY)",
    y = "Precipitation",
    color = "Quantile"
  ) +
  theme_minimal()







