# Generate a stacked bar chart of test history from tests/test_history.csv.
# Run from project root: source("tests/plot_test_history.R")

library(ggplot2)

log_path  <- "tests/test_history.csv"
out_path  <- "tests/test_history.png"

history <- read.csv(log_path, stringsAsFactors = FALSE)
history$run <- seq_len(nrow(history))

# Reshape to long form
long <- rbind(
  data.frame(run = history$run, label = history$label,
             status = "Pass",  n = history$pass),
  data.frame(run = history$run, label = history$label,
             status = "Fail",  n = history$fail),
  data.frame(run = history$run, label = history$label,
             status = "Skip",  n = history$skip),
  data.frame(run = history$run, label = history$label,
             status = "Warn",  n = history$warn)
)
long$status <- factor(long$status, levels = c("Fail", "Warn", "Skip", "Pass"))
long$label  <- factor(long$label, levels = history$label)

colours <- c(Pass = "#2ecc71", Skip = "#f39c12", Warn = "#e67e22", Fail = "#e74c3c")

p <- ggplot(long, aes(x = label, y = n, fill = status)) +
  geom_col(width = 0.65) +
  geom_text(
    data = subset(long, n > 0),
    aes(label = n),
    position = position_stack(vjust = 0.5),
    colour = "white", fontface = "bold", size = 3.5
  ) +
  scale_fill_manual(values = colours, name = "Result") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  labs(
    title    = "Test suite history — SizeShiftsImplicationsV2",
    subtitle = "Each bar = one devtools::test() run",
    x        = NULL,
    y        = "Number of tests"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    legend.position = "top",
    plot.title      = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

ggsave(out_path, p, width = 8, height = 5, dpi = 150)
cat("Saved:", out_path, "\n")
