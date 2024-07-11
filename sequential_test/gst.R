library("ldbounds")
# library("rpact")

########################################
# 1件のGSTシミュレーションを実行する関数
########################################
measure_result <- function(n_max, k_max, pA, pB, critical_values) {
  # 1回のシミュレーション結果の初期化
  result <- list(
    k_max           = k_max,
    k               = k_max,
    t               = c(),
    n_total         = c(),
    z               = c(),
    critical_values = c(),
    rejected        = c(),
    summary         = list()
  )

  # 中間解析あたりのサンプルサイズ。
  # 等分したいため、round()は使わない。
  # ceiling()ではなくfloor()としたのは、GSTが固定サンプルサイズ検定に対して、検出力の観点で有利にならないようにするため。
  n <- floor(n_max / k_max)

  # 各群の累積のサンプルサイズ
  n_total <- 0

  # A, B群それぞれの標本ベクトル
  A <- c()
  B <- c()

  for (k in 1:k_max) {
    # 中間解析で得たサンプルを蓄積する
    A <- c(A, rbinom(n, 1, pA))
    B <- c(B, rbinom(n, 1, pB))
    n_total <- n_total + n

    # 蓄積されたサンプルで検定統計量を計算
    meanA <- mean(A)
    meanB <- mean(B)

    if (meanA == 1 && meanB == 1) {
      Zk <- 0
    } else if (meanA == 0 && meanB == 0) {
      Zk <- 0
    } else {
      Zk <- (meanA - meanB) / sqrt(meanA * (1 - meanA) / n_total + meanB * (1 - meanB) / n_total)
    }

    result$k               <- k
    result$t               <- c(result$t, k / k_max)
    result$n_total         <- c(result$n_total, n_total)
    result$z               <- c(result$z, Zk)
    result$critical_values <- c(result$critical_values, critical_values[k])

    # 棄却域を超えるかどうか確認
    if (abs(Zk) > critical_values[k]) {
      # 帰無仮説を棄却して、検定を終了する
      result$rejected <- c(result$rejected, "Yes")

      break
    }

    result$rejected <- c(result$rejected, "No")
  }

  result$summary <- list(
    k_last   = result$k,
    n_total  = n_total,
    rejected = tail(result$rejected, 1)
  )

  return(result)
}

# シミュレーション用パラメータ設定
set.seed(123)            # 再現のためのシード設定
simulation_times <- 1000 # シミュレーション回数

# 仮説検定用のパラメータ設定
alpha   <- 0.05          # 期待するtype Iエラー率（FPR）
beta    <- 0.20          # 期待するtype IIエラー率
pB      <- 0.95          # B群の母比率
pA_list <- c(0.96, 0.94, 0.99, 0.90) # 仮説ごとのA群の母比率

# GST特有の仮説検定用のパラメータ設定
sampling_size_multiplier <- c(0.1, 0.5, 1.0, 1.5) # 最大サンプルサイズ。固定サンプルサイズに対する割合ごとにシミュレーション。
k_max <- 20              # 中間解析の回数
rho   <- 3               # Kim & DeMets Power family型のエラー消費関数のパラメータ

# エラー消費関数をごとの棄却域
# ldboundsを使う場合
designs <- list(
  list(
    type = "Pocock-like error spending function",
    critical_values = ldBounds(k_max, alpha = alpha, sides = 1, iuse = 2)$upper.bounds
  ),
  list(
    type = "O'Brien-Fleming-like error spending function",
    critical_values = ldBounds(k_max, alpha = alpha, sides = 1, iuse = 1)$upper.bounds
  ),
  list(
    type = "Kim & DeMets Power family error spending function",
    critical_values = ldBounds(k_max, alpha = alpha, sides = 1, iuse = 3, phi = rho)$upper.bounds
  ),
  list(
    type = "Hwang, Shih & DeCani Gammma family with gamma = +2 error spending function",
    critical_values = ldBounds(k_max, alpha = alpha, sides = 1, iuse = 4, phi = 2)$upper.bounds
  ),
  list(
    type = "Hwang, Shih & DeCani Gammma family with gamma = -2 error spending function",
    critical_values = ldBounds(k_max, alpha = alpha, sides = 1, iuse = 4, phi = -2)$upper.bounds
  )
)
# rpactを使う場合
# designs <- list(
#   list(
#     type = "Pocock-like error spending function",
#     critical_values = as.data.frame(getDesignGroupSequential(kMax = k_max, alpha = alpha, beta = beta, sided = 1, typeOfDesign = "asP"))$criticalValues
#   ),
#   list(
#     type = "O'Brien-Fleming-like error spending function",
#     critical_values = as.data.frame(getDesignGroupSequential(kMax = k_max, alpha = alpha, beta = beta, sided = 1, typeOfDesign = "asOF"))$criticalValues
#   ),
#   list(
#     type = "Kim & DeMets Power family error spending function",
#     critical_values = as.data.frame(getDesignGroupSequential(kMax = k_max, alpha = alpha, beta = beta, sided = 1, typeOfDesign = "asKD", gammaA = rho))$criticalValues
#   ),
#   list(
#     type = "Hwang, Shih & DeCani Gammma family with gamma = +2 error spending function",
#     critical_values = as.data.frame(getDesignGroupSequential(kMax = k_max, alpha = alpha, beta = beta, sided = 1, typeOfDesign = "asHSD", gammaA = 2))$criticalValues
#   ),
#   list(
#     type = "Hwang, Shih & DeCani Gammma family with gamma = -2 error spending function",
#     critical_values = as.data.frame(getDesignGroupSequential(kMax = k_max, alpha = alpha, beta = beta, sided = 1, typeOfDesign = "asHSD", gammaA = -2))$criticalValues
#   )
# )

# エラー消費関数ごとにシミュレーションする
for (design in designs) {
  print(design$type)

  # 結果リスト
  summary <- list()

  # 仮説ごとにシミュレーションを実行
  for (pA in pA_list) {
    fixed_sample_size <- power.prop.test(
      p1 = pA,
      p2 = pB,
      sig.level = alpha,
      power = 1 - beta,
      alternative = "one.sided"
    )$n

    # サンプリングサイズの割合ごとにシミュレーションを実行
    for (multiplier in sampling_size_multiplier) {
      k_last_results <- c()
      n_total_results <- c()
      rejected_results <- c()

      n_max <- ceiling(fixed_sample_size * multiplier)

      for (i in 1:simulation_times) {
        result <- measure_result(n_max = n_max, k_max = k_max, pA = pA, pB = pB, critical_values = design$critical_values)
        k_last_results <- c(k_last_results, result$summary$k_last)
        n_total_results <- c(n_total_results, result$summary$n_total)
        rejected_results <- c(rejected_results, result$summary$rejected)
      }

      # 帰無仮説を棄却できた回数
      rejected_times <- sum(rejected_results == "Yes")

      # 平均サンプルサイズ
      n_average <- round(mean(n_total_results))

      summary <- append(summary, list(list(
        pA                       = pA,
        simulation_times         = simulation_times,
        k_average                = mean(k_last_results),# 平均中間解析回数
        sampling_size_multiplier = multiplier,
        n_max                    = n_max,
        n_average                = n_average,
        sample_size_reduction    = (fixed_sample_size - n_average) / fixed_sample_size, # 固定サンプルサイズに対するサンプルサイズの削減率
        rejected_times           = rejected_times,
        power                    = rejected_times / simulation_times # 実際の検出力
      )))
    }
  }

  # シミュレーション結果をテーブルに変換して表示
  df <- as.data.frame(t(sapply(summary, function(x) unlist(x))))
  print(df)
}
