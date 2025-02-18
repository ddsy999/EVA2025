library(dplyr)
library(ggplot2)
library(zoo)

load("data//run1.RData")
load("data//run2.RData")
load("data//run3.RData")
load("data//run4.RData")


# 1. 날짜 벡터 생성 (1850-01-01 ~ 2014-12-31, 총 60225일)
dates <- seq(from = as.Date("1850-01-01"), by = "day", length.out = 60225)

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
