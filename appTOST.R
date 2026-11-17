library(shiny)
library(bslib)
library(ggplot2)

# Préréglages des 4 scénarios du poster (Figure 3)
presets <- list(
  A = list(diff = 0.02,  delta = 0.40, sd = 0.25, n = 30),
  B = list(diff = -0.45, delta = 0.40, sd = 0.30, n = 25),
  C = list(diff = 0.22,  delta = 0.40, sd = 0.18, n = 40),
  D = list(diff = 0.12,  delta = 0.40, sd = 0.55, n = 10)
)

ui <- page_sidebar(
  title = span("Démonstrateur TOST — Master 2 CACQ-OPEx (UBO)", class = "fw-bold"),
  theme = bs_theme(
    version = 5, bootswatch = "flatly",
    primary = "#1e6b43", bg = "#f4f8f5", fg = "#1c2833",
    base_font = font_google("Inter"), heading_font = font_google("Montserrat")
  ),
  
  sidebar = sidebar(
    width = 330, bg = "#ffffff",
    h6("Scénarios du Poster (Figure 3)", class = "fw-bold text-uppercase text-muted mb-2"),
    div(
      class = "d-grid gap-2 mb-3",
      actionButton("btn_A", "Cas A : Équivalent & Non différent", class = "btn-success btn-sm text-start"),
      actionButton("btn_B", "Cas B : Non équivalent & Différent", class = "btn-danger btn-sm text-start"),
      actionButton("btn_C", "Cas C : Équivalent & Différent", class = "btn-info btn-sm text-white text-start"),
      actionButton("btn_D", "Cas D : Indéterminé (Sous-puissant)", class = "btn-warning btn-sm text-white text-start")
    ),
    hr(class = "my-2"),
    h6("Paramètres Expérimentaux", class = "fw-bold text-uppercase text-muted mb-2"),
    sliderInput("diff_mean", "Différence observée (X̄₁ - X̄₂) :", min = -1, max = 1, value = 0.02, step = 0.01),
    sliderInput("delta", "Marge de tolérance (±Δ) :", min = 0.1, max = 1, value = 0.40, step = 0.05),
    numericInput("sd_pooled", "Écart-type combiné (s) :", value = 0.25, min = 0.01, max = 5, step = 0.05),
    numericInput("n_group", "Taille d'échantillon / groupe (n) :", value = 30, min = 3, max = 500, step = 1),
    sliderInput("alpha", "Seuil unilatéral (α) :", min = 0.01, max = 0.10, value = 0.05, step = 0.01)
  ),
  
  div(
    uiOutput("verdict_banner"),
    card(
      card_header(
        class = "fw-bold d-flex justify-content-between align-items-center",
        span("Visualisation Géométrique"),
        span(class = "badge bg-success", "Superposition des IC & Zone d'Équivalence")
      ),
      card_body(plotOutput("plot_intervals", height = "340px"))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Test d'Équivalence (TOST)", class = "fw-bold text-success"), card_body(tableOutput("table_tost"), class = "p-2")),
      card(card_header("Test Bilatéral Classique (Student)", class = "fw-bold text-primary"), card_body(tableOutput("table_nhst"), class = "p-2"))
    )
  )
)

server <- function(input, output, session) {
  
  # Application factorisée des 4 scénarios
  lapply(names(presets), function(id) {
    observeEvent(input[[paste0("btn_", id)]], {
      p <- presets[[id]]
      updateSliderInput(session, "diff_mean", value = p$diff)
      updateSliderInput(session, "delta", value = p$delta)
      updateNumericInput(session, "sd_pooled", value = p$sd)
      updateNumericInput(session, "n_group", value = p$n)
      updateSliderInput(session, "alpha", value = 0.05)
    })
  })
  
  stats_data <- reactive({
    req(input$sd_pooled, input$n_group)
    s <- max(0.001, input$sd_pooled)
    n <- max(3, round(input$n_group))
    diff <- input$diff_mean; delta <- input$delta; alpha <- input$alpha
    
    df <- 2 * n - 2
    se <- s * sqrt(2 / n)
    pct_tost <- round((1 - 2 * alpha) * 100, 1)
    pct_nhst <- round((1 - alpha) * 100, 1)
    
    # Bilatéral classique (NHST)
    t_nhst <- diff / se
    p_nhst <- 2 * (1 - pt(abs(t_nhst), df))
    ci_nhst <- diff + c(-1, 1) * qt(1 - alpha / 2, df) * se
    diff_sig <- p_nhst < alpha
    
    # TOST unilatéral & Intervalle de Westlake
    t1 <- (diff + delta) / se; p1 <- 1 - pt(t1, df)
    t2 <- (diff - delta) / se; p2 <- pt(t2, df)
    crit_tost <- qt(1 - alpha, df)
    ci_tost <- diff + c(-1, 1) * crit_tost * se
    equiv <- (ci_tost[1] > -delta) && (ci_tost[2] < delta)
    power <- max(0, 1 - 2 * pt(crit_tost, df, ncp = delta / se))
    
    title <- if (equiv && !diff_sig) "Cas A : Équivalent & Non différent"
    else if (!equiv && diff_sig) "Cas B : Non équivalent & Différent"
    else if (equiv && diff_sig) "Cas C : Équivalent & Différent (Écart statistique mais négligeable en pratique)"
    else "Cas D : Non équivalent & Non différent (Zone d'incertitude / Manque de puissance)"
    
    list(diff = diff, delta = delta, df = df, alpha = alpha, pct_tost = pct_tost, pct_nhst = pct_nhst,
         t_nhst = t_nhst, p_nhst = p_nhst, ci_nhst = ci_nhst, diff_sig = diff_sig,
         t1 = t1, p1 = p1, t2 = t2, p2 = p2, p_max = max(p1, p2), ci_tost = ci_tost,
         equiv = equiv, power = power, title = title)
  })
  
  output$verdict_banner <- renderUI({
    r <- stats_data()
    bg <- if (r$equiv) "#d4efdf" else "#fadbd8"
    col <- if (r$equiv) "#27ae60" else "#c0392b"
    msg <- if (r$equiv) {
      sprintf("Équivalence validée (max(p1, p2) = %.4f ≤ %.2f). L'IC à %s%% est inclus dans [-Δ ; +Δ].", r$p_max, r$alpha, r$pct_tost)
    } else {
      sprintf("Équivalence rejetée (max(p1, p2) = %.4f > %.2f). L'IC à %s%% dépasse au moins une borne.", r$p_max, r$alpha, r$pct_tost)
    }
    div(
      style = sprintf("background-color:%s; border-left:6px solid %s; padding:12px 16px; border-radius:6px; margin-bottom:14px;", bg, col),
      h5(r$title, style = sprintf("color:%s; font-weight:700; margin-bottom:4px;", col)),
      div(class = "small", tags$b("Décision TOST : "), msg, tags$span(" | "), tags$b("Puissance théorique (à diff=0) : "), sprintf("%.1f%%", r$power * 100))
    )
  })
  
  output$plot_intervals <- renderPlot({
    r <- stats_data()
    lbl <- c(sprintf("IC %s%% (Student)", r$pct_nhst), sprintf("IC %s%% (TOST)", r$pct_tost))
    df_p <- data.frame(
      y = c(0.7, 1.3), x = r$diff,
      xmin = c(r$ci_nhst[1], r$ci_tost[1]), xmax = c(r$ci_nhst[2], r$ci_tost[2]),
      col = c(ifelse(r$diff_sig, "#2980b9", "#7f8c8d"), ifelse(r$equiv, "#27ae60", "#c0392b"))
    )
    xlim <- max(abs(c(r$diff, r$delta, r$ci_tost, r$ci_nhst))) * 1.3
    
    ggplot() +
      annotate("rect", xmin = -r$delta, xmax = r$delta, ymin = 0.1, ymax = 1.9, fill = "#27ae60", alpha = 0.15) +
      geom_vline(xintercept = c(-r$delta, r$delta), linetype = "dashed", color = "#1e6b43", linewidth = 1) +
      geom_vline(xintercept = 0, color = "#2c3e50", linewidth = 0.8) +
      annotate("text", x = c(-r$delta, r$delta), y = 1.82, label = c("-Δ", "+Δ"), color = "#1e6b43", fontface = "bold", hjust = c(1.3, -0.3), size = 4.5) +
      geom_errorbar(data = df_p, aes(y = y, xmin = xmin, xmax = xmax, color = col), width = 0.18, linewidth = 1.6, orientation = "y") +
      geom_point(data = df_p, aes(x = x, y = y, color = col), size = 4.5, shape = 21, fill = "white", stroke = 2.4) +
      scale_color_identity() +
      scale_y_continuous(breaks = c(0.7, 1.3), labels = lbl, limits = c(0.15, 1.9)) +
      coord_cartesian(xlim = c(-xlim, xlim)) +
      labs(title = "Comparaison géométrique des Intervalles de Confiance",
           subtitle = sprintf("L'équivalence requiert l'inclusion stricte de l'IC à %s%% dans [-Δ ; +Δ]", r$pct_tost),
           x = "Différence de moyenne (X̄₁ - X̄₂)", y = NULL) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), axis.text.y = element_text(face = "bold", size = 12),
            plot.title = element_text(face = "bold", size = 14, color = "#1e6b43"), plot.subtitle = element_text(color = "#566573", size = 11))
  })
  
  output$table_tost <- renderTable({
    r <- stats_data()
    data.frame(
      Paramètre = c("Test bas t1 (-Δ)", "Test haut t2 (+Δ)", "Critère max(p1, p2)", sprintf("IC %s%% (TOST)", r$pct_tost)),
      Résultat = c(sprintf("t1 = %.2f (p1 = %.4f)", r$t1, r$p1), sprintf("t2 = %.2f (p2 = %.4f)", r$t2, r$p2),
                   sprintf("%.4f (seuil α = %.2f)", r$p_max, r$alpha), sprintf("[ %.3f ; %.3f ]", r$ci_tost[1], r$ci_tost[2]))
    )
  }, striped = TRUE, hover = TRUE, width = "100%")
  
  output$table_nhst <- renderTable({
    r <- stats_data()
    data.frame(
      Paramètre = c("Statistique t", "Degrés de liberté (ddl)", "p-value bilatérale", sprintf("IC %s%% (Student)", r$pct_nhst)),
      Résultat = c(sprintf("t = %.2f", r$t_nhst), sprintf("%d", r$df), sprintf("%.4f (seuil α = %.2f)", r$p_nhst, r$alpha),
                   sprintf("[ %.3f ; %.3f ]", r$ci_nhst[1], r$ci_nhst[2]))
    )
  }, striped = TRUE, hover = TRUE, width = "100%")
}

shinyApp(ui = ui, server = server)