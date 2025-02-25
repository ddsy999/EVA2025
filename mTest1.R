

################# X_t , X_{t-1} plot ##############
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







################# X_t , X_{t-1} auto  ##############
# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run1[1,1,]

# (1) ACF를 통한 자기상관 확인
acf(data_series, lag.max=10, main="ACF: 시차 1~10") 

# (2) Ljung-Box 검정 (lag=1)
Box.test(data_series, lag=1, type="Ljung-Box")

# (3) 선형 회귀 검정 (X_{t-1} → X_t)
df <- data.frame(X_t = data_series[2:length(data_series)],
                 X_t1 = data_series[1:(length(data_series)-1)])

lm_result <- lm(X_t ~ X_t1, data=df)
summary(lm_result)












################# X_{t-1} 의 크기에 따른 평균 변화량 ##############
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




################# X_{t-1} 의 크기에 따른 계수 #############
# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run1[1,1,]

# X_{t-1}과 X_t 정의 (로그 변환 없이)
X_t1 <- data_series[1:(length(data_series)-1)]
X_t <- data_series[2:length(data_series)]

# 데이터프레임 생성
df <- data.frame(X_t1, X_t)

# X_{t-1}을 분위수(quantile) 기준으로 10개 그룹으로 나눔
df <- df %>%
  mutate(X_bin = ntile(X_t1, 10))  # 10분위로 구분

# 각 구간별로 AR(1) 모델 적합 및 계수, p-value 계산
ar_results <- df %>%
  group_by(X_bin) %>%
  summarise(beta1 = coef(lm(X_t ~ X_t1, data = cur_data()))[2],
            p_value = summary(lm(X_t ~ X_t1, data = cur_data()))$coefficients[2,4])

# 결과 출력
print(ar_results)

# 계수 시각화
ggplot(ar_results, aes(x=factor(X_bin), y=beta1)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  labs(title="X_{t-1} 크기에 따른 AR(1) 계수 변화 (원 자료)", 
       x="X_{t-1} 구간", y="AR(1) 계수") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






################# X_{t-1} 의 크기에 따른 X_t 값 산점도  #############
# 특정 위치 (1,1)에서 시계열 데이터 추출
data_series <- run3[1,1,]

# X_{t-1}과 X_t 생성
X_t1 <- data_series[1:(length(data_series)-1)]
X_t <- data_series[2:length(data_series)]

# 데이터프레임 생성
df <- data.frame(X_t1, X_t)

# X_{t-1}을 분위수(quantile) 기준으로 10개 그룹으로 나눔
df <- df %>%
  mutate(X_bin = ntile(X_t1, 10))  # 10분위(Decile)로 구분

# facet_wrap을 활용한 그룹별 산점도
ggplot(df, aes(x = X_t1, y = X_t)) +
  geom_point(alpha = 0.6, color = "blue") +  # 반투명한 파란색 점
  theme_minimal() +
  labs(title = "구간별 X_{t-1} vs. X_t 산점도 (Facet Free)",
       x = "X_{t-1} (이전 시점 값)",
       y = "X_t (현재 시점 값)") +
  facet_wrap(~X_bin, scales = "free")  # 그룹별로 분리 및 자유 스케일 적용



