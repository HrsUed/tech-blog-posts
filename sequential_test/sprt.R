######################################
# 1件のシミュレーションを実行する関数
######################################
measure_result <- function(alpha, beta, p0, p1) {
  q0 <- 1 - p0
  q1 <- 1 - p1
  c1 <- beta / (1 - alpha)
  c2 <- (1 - beta) / alpha

  # 固定サンプルサイズの見積もり
  # See https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/power.prop.test
  fixed_size <- power.prop.test(
    p1 = p0,
    p2 = p1,
    sig.level = alpha,
    power = (1-beta),
    alternative = "one.sided"
  )$n 
  
  # 対数尤度関数の初期値
  L <- 0
  
  # 結果の初期化
  result <- list(
    null = FALSE,
    alternative = FALSE,
    size = fixed_size,
    fixed_size = fixed_size
  )
  
  # SPRTによる逐次検定
  # 終わらない可能性を考慮して、固定サンプルサイズを上限とする
  for (i in 1:fixed_size) {
    # ベルヌーイ試行
    x <- rbinom(1, 1, p1)
    
    # 尤度関数の更新
    L <- L + log(q1/q0) + x * log(p1*q0/p0/q1)

    # 検定結果の確認
    if (L < log(c1) || log(c2) < L) {
      result$null        <- (L < log(c1))
      result$alternative <- (L > log(c2))
      result$size        <- i

      break
    }
  }
  
  return(result)
}

# シミュレーション用パラメータ設定
sim_times <- 10000 # シミュレーション回数
set.seed(123)      # 再現性のためのシード設定
p0    <- 0.95      # 帰無仮説
alpha <- 0.05      # 期待するtype Iエラー率（FPR）
beta  <- 0.20      # 期待するtype IIエラー率
p1_list <- c(0.96, 0.94, 0.99, 0.90) # シミュレーションする対立仮説の値

# 結果リスト
summary <- list()

# 対立仮説ごとにシミュレーションを実行
for (p1 in p1_list) {
  results <- list()
  
  for (i in 1:sim_times) {
    result <- measure_result(alpha = alpha, beta = beta, p0 = p0, p1 = p1)
    results <- append(results, list(result))
  }
  
  summary <- append(
    summary,
    list(list(
      alternative = p1,
      sim_times = sim_times,
      mean = round(mean(sapply(results, function(x) { x$size }))),
      sd = round(sd(sapply(results, function(x) { x$size }))),
      power = sum(sapply(results, function(x) { x$alternative })) / sim_times,
      stop = sum(sapply(results, function(x) { !x$null & !x$alternative })) / sim_times
    ))
  )
}

# シミュレーション結果をテーブルに変換
df <- as.data.frame(t(sapply(summary, function(x) unlist(x))))
print(df)