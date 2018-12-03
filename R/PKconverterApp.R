#' Shiny App for converting the pharmacokinetic parameters
#'
#' Pharmacokinetics is the study of drug absorption, distribution, metabolism,
#' and excretion. The pharmacokinetics model explains that how the drug
#' concentration change as the drug moves through the different compartments
#' of the body. For pharmacokinetic modeling and analysis, it is essential
#' to understand the basic pharmacokinetic parameters. All parameters are
#' considered, but only some of parameters are used in the model.
#' Therefore,we need to convert the estimated parameters to the other
#' parameters after fitting the specific pharmacokinetic model.
#' For Shiny App on the web, visit
#' \href{https://ek-lee.shinyapps.io/PKconverter/}{PKconverter shiny app},
#' @usage PKconverterApp()
#' @examples
#' #PKconverterApp()
#' @export
#'
PKconverterApp<-function(){

  PKconverter.env<-new.env()
  PK_app=shiny::shinyApp(
    ui=shiny::fluidPage(
      theme=shinythemes::shinytheme("cerulean"),
      shiny::shinyUI(
        shiny::navbarPage("PK Parameter Converter",

        ############################sheet1############################

        shiny::tabPanel(title="Model 1",
          shiny::sidebarLayout(
           shiny::column(width=5,
            shinydashboard::box(title="Select your model",width=NULL,
              solidHeader=TRUE,status="warning",
              background="light-blue",
              shiny::radioButtons("model1", "MODEL TYPE :",
                list("One compartment model" ="one",
                  "Two compartment model" ="two",
                  "Three compartment model" ="three"))),
            shiny::tags$hr(),
            shiny::tags$hr(),
            shinydashboard::box(title="Enter your estimate and std.err",
              width=NULL,solidHeader=TRUE,status="warning",
              background="light-blue",
              shiny::tabsetPanel(id="input1",type="pills",
                shiny::tabPanel("Estimate / Std.err", value="meanandsd",
                  shiny::tags$hr(),
                  shiny::conditionalPanel(condition="input.model1=='one'",
                    shiny::uiOutput("comp1_1_1"),
                    shiny::uiOutput("comp1_1_2")),
                  shiny::conditionalPanel(condition="input.model1=='two'",
                    shiny::uiOutput("comp2_1_1"),
                    shiny::uiOutput("comp2_1_2"),
                    shiny::uiOutput("comp2_1_3"),
                    shiny::uiOutput("comp2_1_4")),
                  shiny::conditionalPanel(condition="input.model1=='three' ",
                    shiny::uiOutput("comp3_1_1"),
                    shiny::uiOutput("comp3_1_2"),
                    shiny::uiOutput("comp3_1_3"),
                    shiny::uiOutput("comp3_1_4"),
                    shiny::uiOutput("comp3_1_5"),
                    shiny::uiOutput("comp3_1_6"))),

                shiny::tabPanel("Covariance", value="varcov",
                  shiny::conditionalPanel(condition="input.model1=='one'||
                    input.model1=='two'||input.model1=='three'",
                    shiny::fluidRow(shiny::column(4,
                      shiny::uiOutput("v1cl1_1")))),
                  shiny::conditionalPanel(condition="input.model1=='two'||
                    input.model1=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v2cl1_1")),
                    shiny::column(4,shiny::uiOutput("v1v2_1")))),
                  shiny::conditionalPanel(condition="input.model1=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v3cl1_1")),
                    shiny::column(4,shiny::uiOutput("v1v3_1")),
                    shiny::column(4,shiny::uiOutput("v2v3_1")))),
                  shiny::conditionalPanel(condition="input.model1=='two'||
                    input.model1=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v1cl2_1")),
                    shiny::column(4,shiny::uiOutput("v2cl2_1")),
                    shiny::column(4,shiny::uiOutput("cl1cl2_1")))),
                  shiny::conditionalPanel(condition="input.model1=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v1cl3_1")),
                    shiny::column(4,shiny::uiOutput("v3cl2_1")),
                    shiny::column(4,shiny::uiOutput("cl1cl3_1")))),
                  shiny::conditionalPanel(condition="input.model1=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v2cl3_1")),
                    shiny::column(4,shiny::uiOutput("v3cl3_1")),
                    shiny::column(4,shiny::uiOutput("cl2cl3_1"))))))),
            shiny::tags$hr(),
            shiny::tags$hr(),
            shinydashboard::box(title="Save results as a file",width=NULL,
              solidHeader=TRUE,status="warning",
              background="light-blue",
              shiny::radioButtons("type1", "File type:",
                choices = c("Excel (CSV)" ,
                  "Text (tab separated)" ,
                  "Text (Space Separated)")),
              shiny::conditionalPanel(condition = "input.model1=='one'",
                shiny::downloadButton("downloadData1_1",
                  "Save results to file")),
              shiny::conditionalPanel(condition = "input.model1=='two'",
                shiny::downloadButton("downloadData1_2",
                  "Save results to file")),
              shiny::conditionalPanel(condition = "input.model1=='three'",
                shiny::downloadButton("downloadData1_3",
                  "Save results to file")))),
           shiny::column(width=7,
            shiny::wellPanel(
              shiny::h2("Model 1: Volumes and Clearances"),
                shiny::tags$br(),
                shiny::conditionalPanel(condition = "input.model1=='one'",
                  shiny::h3("One compartment model"),
                  shiny::tableOutput("T1_1")),
                shiny::conditionalPanel(condition = "input.model1=='two'",
                  shiny::h3("Two compartment model"),
                  shiny::tableOutput("T1_2")),
                shiny::conditionalPanel(condition = "input.model1=='three'",
                  shiny::h3("Three compartment model"),
                  shiny::tableOutput("T1_3")))))),

        ############################sheet2############################

        shiny::tabPanel(title="Model 2",
          shiny::sidebarLayout(
           shiny::column(width=5,
            shinydashboard::box(title="Select your model",width=NULL,
              solidHeader=TRUE,status="warning",
              background="light-blue",
              shiny::radioButtons("model2", "MODEL TYPE :",
                list("One compartment model" ="one",
                  "Two compartment model" ="two",
                  "Three compartment model" ="three"))),
            shiny::tags$hr(),
            shiny::tags$hr(),
            shinydashboard::box(title="Enter your estimate and std.err",
              width=NULL,solidHeader=TRUE,status="warning",
              background="light-blue",
              shiny::tabsetPanel(id="input2",
                shiny::tabPanel("Estimate / Std.err", value="meanandsd2",
                  shiny::tags$hr(),
                  shiny::conditionalPanel(condition="input.model2=='one'",
                    shiny::uiOutput("comp1_2_1"),
                    shiny::uiOutput("comp1_2_2")),
                  shiny::conditionalPanel(condition="input.model2=='two'",
                    shiny::uiOutput("comp2_2_1"),
                    shiny::uiOutput("comp2_2_2"),
                    shiny::uiOutput("comp2_2_3"),
                    shiny::uiOutput("comp2_2_4")),
                  shiny::conditionalPanel(condition="input.model2=='three'",
                    shiny::uiOutput("comp3_2_1"),
                    shiny::uiOutput("comp3_2_2"),
                    shiny::uiOutput("comp3_2_3"),
                    shiny::uiOutput("comp3_2_4"),
                    shiny::uiOutput("comp3_2_5"),
                    shiny::uiOutput("comp3_2_6"))),

                shiny::tabPanel("Covariance", value="varcov2",
                  shiny::conditionalPanel(
                    condition="input.model2=='one'||input.model2=='two'
                    ||input.model2=='three'",
                    shiny::fluidRow(shiny::column(4,
                      shiny::uiOutput("v1k10_2")))),
                  shiny::conditionalPanel(
                    condition="input.model2=='two'||input.model2=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v1k12_2")),
                      shiny::column(4,shiny::uiOutput("v1k21_2")))),
                  shiny::conditionalPanel(condition="input.model2=='three'",
                    shiny::fluidRow(shiny::column(4,shiny::uiOutput("v1k13_2")),
                      shiny::column(4,shiny::uiOutput("v1k31_2")),
                      shiny::column(4,shiny::uiOutput("k13k31_2")))),
                  shiny::conditionalPanel(
                    condition="input.model2=='two'||input.model2=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("k10k12_2")),
                        shiny::column(4,shiny::uiOutput("k10k21_2")),
                        shiny::column(4,shiny::uiOutput("k12k21_2")))),
                  shiny::conditionalPanel(condition="input.model2=='three'",
                    shiny::fluidRow(shiny::column(4,
                      shiny::uiOutput("k10k13_2")),
                      shiny::column(4,shiny::uiOutput("k10k31_2")),
                      shiny::column(4,shiny::uiOutput("k12k31_2"))),
                  shiny::conditionalPanel(condition="input.model2=='three'",
                    shiny::fluidRow(shiny::column(4,
                      shiny::uiOutput("k21k13_2")),
                      shiny::column(4,shiny::uiOutput("k21k31_2")),
                      shiny::column(4,shiny::uiOutput("k12k13_2")))))))),
            shiny::tags$hr(),
            shiny::tags$hr(),
            shinydashboard::box(title="Save the result as a file",width=NULL,
              solidHeader=TRUE,status="warning",
              background="light-blue",
              shiny::radioButtons("type2", "File type:",
                choices = c("Excel (CSV)" ,
                  "Text (tab separated)" ,
                  "Text (Space Separated)")),

              shiny::conditionalPanel(condition = "input.model2=='one'",
                shiny::downloadButton("downloadData2_1",
                  "Save results to file")),
              shiny::conditionalPanel(condition = "input.model2=='two'",
                shiny::downloadButton("downloadData2_2",
                  "Save results to file")),
              shiny::conditionalPanel(condition = "input.model2=='three'",
                shiny::downloadButton("downloadData2_3",
                  "Save results to file")))),
           shiny::column(width=7,
            shiny::wellPanel(
              shiny::h2("Model 2: V1, Rate Constants"),
              shiny::tags$br(),
              shiny::conditionalPanel(condition = "input.model2=='one'",
                shiny::h3("One compartment model"),
                shiny::tableOutput("T2_1")),
              shiny::conditionalPanel(condition = "input.model2=='two'",
                shiny::h3("Two compartment model"),
                shiny::tableOutput("T2_2")),
              shiny::conditionalPanel(condition = "input.model2=='three'",
                shiny::h3("Three compartment model"),
                shiny::tableOutput("T2_3")))))),

        ############################sheet3############################

        shiny::tabPanel(title="Model 3",
          shiny::sidebarLayout(
            shiny::column(width=5,
              shinydashboard::box(title="Select your model",width=NULL,
                solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("model3", "MODEL TYPE :",
                  list("One compartment model" ="one",
                    "Two compartment model" ="two",
                    "Three compartment model" ="three"))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Enter your estimate and std.err",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::tabsetPanel(id="input3",
                  shiny::tabPanel("Estimate / Std.err", value="meanandsd3",
                    shiny::conditionalPanel(condition="input.model3=='one' ",
                      shiny::uiOutput("comp1_3_1"),
                      shiny::uiOutput("comp1_3_2")),
                    shiny::conditionalPanel(condition="input.model3=='two'",
                      shiny::uiOutput("comp2_3_1"),
                      shiny::uiOutput("comp2_3_2"),
                      shiny::uiOutput("comp2_3_3"),
                      shiny::uiOutput("comp2_3_4")),
                    shiny::conditionalPanel(condition="input.model3=='three' ",
                      shiny::uiOutput("comp3_3_1"),
                      shiny::uiOutput("comp3_3_2"),
                      shiny::uiOutput("comp3_3_3"),
                      shiny::uiOutput("comp3_3_4"),
                      shiny::uiOutput("comp3_3_5"),
                      shiny::uiOutput("comp3_3_6"))),
                  shiny::tabPanel("Covariance", value="varcov3",
                    shiny::conditionalPanel(condition="input.model3=='one'||
                      input.model3=='two'||input.model3=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("cl1talpha_3")))),
                    shiny::conditionalPanel(condition="input.model3=='two'||
                      input.model3=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("v1talpha_3")),
                        shiny::column(4,shiny::uiOutput("v1tbeta_3")))),
                    shiny::conditionalPanel(condition="input.model3=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("vdtalpha_3")),
                        shiny::column(4,shiny::uiOutput("vdtbeta_3")),
                        shiny::column(4,shiny::uiOutput("vdtgamma_3")))),
                    shiny::conditionalPanel(condition="input.model3=='two'||
                      input.model3=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("v1cl1_3")),
                        shiny::column(4,shiny::uiOutput("cl1tbeta_3")),
                        shiny::column(4,shiny::uiOutput("talphatbeta_3")))),
                    shiny::conditionalPanel(condition="input.model3=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("v1vd_3")),
                        shiny::column(4,shiny::uiOutput("v1tgamma_3")),
                        shiny::column(4,shiny::uiOutput("vdcl1_3")))),
                    shiny::conditionalPanel(condition="input.model3=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("cl1tgamma_3")),
                        shiny::column(4,shiny::uiOutput("talphatgamma_3")),
                        shiny::column(4,shiny::uiOutput("tbetatgamma_3"))))))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Save the result as a file",width=NULL,
                solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("type3", "File type:",
                  choices = c("Excel (CSV)" ,
                    "Text (tab separated)" ,
                    "Text (Space Separated)")),

                shiny::conditionalPanel(condition = "input.model3=='one'",
                  shiny::downloadButton("downloadData3_1",
                    "Save results to file")),
                shiny::conditionalPanel(condition = "input.model3=='two'",
                  shiny::downloadButton("downloadData3_2",
                    "Save results to file")),
                shiny::conditionalPanel(condition = "input.model3=='three'",
                  shiny::downloadButton("downloadData3_3",
                    "Save results to file")))),

            shiny::column(width=7,
              shiny::wellPanel(
                shiny::h2("Model 3: V1, Vdss, Cl, half-lives"),
                  shiny::tags$br(),
                  shiny::conditionalPanel(condition = "input.model3=='one'",
                    shiny::h3("One compartment model"),
                    shiny::tableOutput("T3_1")),
                  shiny::conditionalPanel(condition = "input.model3=='two'",
                    shiny::h3("Two compartment model"),
                    shiny::tableOutput("T3_2")),
                  shiny::conditionalPanel(condition = "input.model3=='three'",
                    shiny::h3("Three compartment model"),
                    shiny::tableOutput("T3_3")))))),

        ############################sheet4############################

        shiny::tabPanel(title="Model 4",
          shiny::sidebarLayout(
            shiny::column(width=5,
              shinydashboard::box(title="Select your model",width=NULL,
                solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("model4", "MODEL TYPE :",
                  list("One compartment model" ="one",
                    "Two compartment model" ="two",
                    "Three compartment model" ="three"))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Enter your estimate and std.err",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::tabsetPanel(id="input4",
                  shiny::tabPanel("Estimate / Std.err", value="meanandsd4",
                    shiny::conditionalPanel(condition="input.model4=='one' ",
                      shiny::uiOutput("comp1_4_1"),
                      shiny::uiOutput("comp1_4_2")),
                    shiny::conditionalPanel(condition="input.model4=='two' ",
                      shiny::uiOutput("comp2_4_1"),
                      shiny::uiOutput("comp2_4_2"),
                      shiny::uiOutput("comp2_4_3"),
                      shiny::uiOutput("comp2_4_4")),
                    shiny::conditionalPanel(condition="input.model4=='three' ",
                      shiny::uiOutput("comp3_4_1"),
                      shiny::uiOutput("comp3_4_2"),
                      shiny::uiOutput("comp3_4_3"),
                      shiny::uiOutput("comp3_4_4"),
                      shiny::uiOutput("comp3_4_5"),
                      shiny::uiOutput("comp3_4_6"))),
                  shiny::tabPanel("Covariance", value="varcov4",
                    shiny::conditionalPanel(condition="input.model4=='one'||
                      input.model4=='two'||input.model4=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("aalpha_4")))),
                    shiny::conditionalPanel(condition="input.model4=='two'||
                      input.model4=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("abeta_4")),
                        shiny::column(4,shiny::uiOutput("ab_4")))),
                    shiny::conditionalPanel(condition="input.model4=='two'||
                      input.model4=='three'",
                        shiny::fluidRow(shiny::column(4,
                          shiny::uiOutput("balpha_4")),
                          shiny::column(4,shiny::uiOutput("bbeta_4")),
                          shiny::column(4,shiny::uiOutput("alphabeta_4")))),
                    shiny::conditionalPanel(condition="input.model4=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("agamma_4")),
                        shiny::column(4,shiny::uiOutput("bgamma_4")),
                        shiny::column(4,shiny::uiOutput("cgamma_4")))),
                    shiny::conditionalPanel(condition="input.model4=='three'",
                      shiny::fluidRow(shiny::column(4,shiny::uiOutput("ac_4")),
                        shiny::column(4,shiny::uiOutput("bc_4")),
                        shiny::column(4,shiny::uiOutput("betagamma_4")))),
                    shiny::conditionalPanel(condition="input.model4=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("alphagamma_4")),
                        shiny::column(4,shiny::uiOutput("calpha_4")),
                        shiny::column(4,shiny::uiOutput("cbeta_4"))))))),
                shiny::tags$hr(),
                shiny::tags$hr(),
                shinydashboard::box(title="Save the result as a file",
                  width=NULL,solidHeader=TRUE,status="warning",
                  background="light-blue",
                  shiny::radioButtons("type4", "File type:",
                    choices = c("Excel (CSV)" ,
                      "Text (tab separated)" ,
                      "Text (Space Separated)")),
                  shiny::conditionalPanel(condition = "input.model4=='one'",
                    shiny::downloadButton("downloadData4_1",
                      "Save results to file")),
                  shiny::conditionalPanel(condition = "input.model4=='two'",
                    shiny::downloadButton("downloadData4_2",
                      "Save results to file")),
                  shiny::conditionalPanel(condition = "input.model4=='three'",
                    shiny::downloadButton("downloadData4_3"
                      ,"Save results to file")))),
            shiny::column(width=7,
              shiny::wellPanel(
                shiny::h2("Model 4: Coefficients and Exponents"),
                shiny::tags$br(),
                shiny::conditionalPanel(condition = "input.model4=='one'",
                  shiny::h3("One compartment model"),
                  shiny::tableOutput("T4_1")),
                shiny::conditionalPanel(condition = "input.model4=='two'",
                  shiny::h3("Two compartment model"),
                  shiny::tableOutput("T4_2")),
                shiny::conditionalPanel(condition = "input.model4=='three'",
                  shiny::h3("Three compartment model"),
                  shiny::tableOutput("T4_3")))))),

        ############################sheet5############################

        shiny::tabPanel(title="Model 5",
          shiny::sidebarLayout(
            shiny::column(width=5,
              shinydashboard::box(title="Select your model",width=NULL,
                solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("model5", "MODEL TYPE :",
                  list("One compartment model" ="one",
                    "Two compartment model" ="two",
                    "Three compartment model" ="three"))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Enter your estimate and std.err",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::tabsetPanel(id="input5",
                  shiny::tabPanel("Estimate / Std.err", value="meanandsd5",
                    shiny::conditionalPanel(condition="input.model5=='one' ",
                      shiny::uiOutput("comp1_5_1"),
                      shiny::uiOutput("comp1_5_2")),
                    shiny::conditionalPanel(condition="input.model5=='two' ",
                      shiny::uiOutput("comp2_5_1"),
                      shiny::uiOutput("comp2_5_2"),
                      shiny::uiOutput("comp2_5_3"),
                      shiny::uiOutput("comp2_5_4")),
                    shiny::conditionalPanel(condition="input.model5=='three' ",
                      shiny::uiOutput("comp3_5_1"),
                      shiny::uiOutput("comp3_5_2"),
                      shiny::uiOutput("comp3_5_3"),
                      shiny::uiOutput("comp3_5_4"),
                      shiny::uiOutput("comp3_5_5"),
                      shiny::uiOutput("comp3_5_6"))),
                  shiny::tabPanel("Covariance", value="varcov5",
                    shiny::conditionalPanel(condition="input.model5=='one'||
                      input.model5=='two'||input.model5=='three'",
                       shiny::fluidRow(shiny::column(4,
                         shiny::uiOutput("v1alpha_5")))),
                    shiny::conditionalPanel(condition="input.model5=='two'||
                      input.model5=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("v1beta_5")),
                        shiny::column(4,shiny::uiOutput("v1k21_5")))),
                    shiny::conditionalPanel(condition="input.model5=='two'||
                      input.model5=='three'",
                        shiny::fluidRow(shiny::column(4,
                          shiny::uiOutput("alphabeta_5")),
                          shiny::column(4,shiny::uiOutput("alphak21_5")),
                          shiny::column(4,shiny::uiOutput("betak21_5")))),
                    shiny::conditionalPanel(condition="input.model5=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("alphak31_5")),
                        shiny::column(4,shiny::uiOutput("betak31_5")),
                        shiny::column(4,shiny::uiOutput("gammak31_5")))),
                    shiny::conditionalPanel(condition="input.model5=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("v1gamma_5")),
                        shiny::column(4,shiny::uiOutput("v1k31_5")),
                        shiny::column(4,shiny::uiOutput("alphagamma_5")))),
                    shiny::conditionalPanel(condition="input.model5=='three'",
                      shiny::fluidRow(shiny::column(4,
                        shiny::uiOutput("betagamma_5")),
                        shiny::column(4,shiny::uiOutput("gammak21_5")),
                        shiny::column(4,shiny::uiOutput("k21k31_5"))))))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Save the result as a file",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("type5", "File type:",
                  choices = c("Excel (CSV)" ,
                    "Text (tab separated)" ,
                    "Text (Space Separated)")),
                shiny::conditionalPanel(condition = "input.model5=='one'",
                  shiny::downloadButton("downloadData5_1",
                    "Save results to file")),
                shiny::conditionalPanel(condition = "input.model5=='two'",
                  shiny::downloadButton("downloadData5_2",
                    "Save results to file")),
                shiny::conditionalPanel(condition = "input.model5=='three'",
                  shiny::downloadButton("downloadData5_3",
                    "Save results to file")))),

              shiny::column(width=7,
                shiny::wellPanel(
                  shiny::h2("Model 5: V1, Exponents, K21, K31"),
                  shiny::tags$br(),
                  shiny::conditionalPanel(condition = "input.model5=='one'",
                    shiny::h3("One compartment model"),
                    shiny::tableOutput("T5_1")),
                  shiny::conditionalPanel(condition = "input.model5=='two'",
                    shiny::h3("Two compartment model"),
                    shiny::tableOutput("T5_2")),
                  shiny::conditionalPanel(condition = "input.model5=='three'",
                    shiny::h3("Three compartment model"),
                    shiny::tableOutput("T5_3")))))),

        ######################### Indiv parameter converter ##################

        shiny::tabPanel("Indiv. Parameter Converter",
          shiny::fluidRow(
            shiny::column(6,
              shinydashboard::box(title="Select your compartment model",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("modelF1", "MODEL TYPE :",
                  list("One compartment model" ="one",
                    "Two compartment model" ="two",
                    "Three compartment model" ="three"))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Select your input data type",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::radioButtons("modelF2", "DATA TYPE :",
                  list("Model 1: Volumes and Clearances" ="M1",
                    "Model 2: V1, Rate Constants" ="M2",
                    "Model 3: V1, Vdss, Cl, half-lives" ="M3",
                    "Model 4: Coefficients and Exponents" ="M4",
                    "Model 5: V1,Exponents,K21,K31" ="M5"))),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Open your individual PK parameters",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::fileInput("inputfile",
                  label="File should be saved as CSV")),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shinydashboard::box(title="Match PK parameters",
                width=NULL,solidHeader=TRUE,status="warning",
                background="light-blue",
                shiny::uiOutput("choose_PK"))),
            shiny::column(6,
              shiny::downloadButton("SaveIndivResult",
                "Save individual parameters to file"),
              shiny::tags$hr(),
              shiny::tags$hr(),
              shiny::tableOutput("ShowResult"))))))), #end ui

    ############################ SERVER  ####################################

    server<-function(input,output,session){

      ## sheet1 ########
      ###### INPUT #####

      output$comp1_1_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean1_1",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd1_1",
            label = "V1 Std.err")))})

      output$comp1_1_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Clearance_mean1_1",
            label = "Cl1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Clearance_sd1_1",
            label = "Cl1 Std.err")))})

      output$comp2_1_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean1_2",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd1_2",
            label = "V1 Std.err")))})

      output$comp2_1_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean2_2",
            label = "V2  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd2_2",
            label = "V2 Std.err")))})

      output$comp2_1_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Clearance_mean1_2",
            label = "Cl1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Clearance_sd1_2",
            label = "Cl1 Std.err")))})

      output$comp2_1_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Clearance_mean2_2",
            label = "Cl2  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Clearance_sd2_2",
            label = "Cl2 Std.err")))})

      output$comp3_1_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean1_3",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd1_3",
            label = "V1 Std.err")))})

      output$comp3_1_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean2_3",
            label = "V2  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd2_3",
            label = "V2 Std.err")))})

      output$comp3_1_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean3_3",
            label = "V3  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd3_3",
            label = "V3 Std.err")))})

      output$comp3_1_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
        shiny::textInput(inputId = "Clearance_mean1_3",
           label = "Cl1  Estimate")),
         shiny::column(6,shiny::textInput(inputId = "Clearance_sd1_3",
           label = "Cl1 Std.err")))})

      output$comp3_1_5<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
         shiny::textInput(inputId = "Clearance_mean2_3",
           label = "Cl2  Estimate")),
         shiny::column(6,shiny::textInput(inputId = "Clearance_sd2_3",
           label = "Cl2 Std.err")))})

      output$comp3_1_6<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Clearance_mean3_3",
            label = "Cl3  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Clearance_sd3_3",
            label = "Cl3 Std.err")))})

      ##### covariance #####
      ##### One compartment #####

      output$v1cl1_1<-shiny::renderUI({
        shiny::textInput(inputId = "v1cl1_1",
                         label = "Cov(V1,Cl1)", value = 0)})

      ##### Two compartment #####

      output$v1cl2_1<-shiny::renderUI({
        shiny::textInput(inputId = "v1cl2_1",
                         label = "Cov(V1,Cl2)", value = 0)})

      output$v2cl1_1<-shiny::renderUI({
        shiny::textInput(inputId = "v2cl1_1",
                         label = "Cov(V2,Cl1)", value = 0)})

      output$v2cl2_1<-shiny::renderUI({
        shiny::textInput(inputId = "v2cl2_1",
                         label = "Cov(V2,Cl2)", value = 0)})

      output$v1v2_1<-shiny::renderUI({
        shiny::textInput(inputId = "v1v2_1",
                         label = "Cov(V1,V2)", value = 0)})

      output$cl1cl2_1<-shiny::renderUI({
        shiny::textInput(inputId = "cl1cl2_1",
                         label = "Cov(Cl1,Cl2)", value = 0)})

      ##### Three compartment #####

      output$v1cl3_1<-shiny::renderUI({
        shiny::textInput(inputId = "v1cl3_1",
                         label = "Cov(V1,Cl3)", value = 0)})

      output$v2cl3_1<-shiny::renderUI({
        shiny::textInput(inputId = "v2cl3_1",
                         label = "Cov(V2,Cl3)", value = 0)})

      output$v3cl1_1<-shiny::renderUI({
        shiny::textInput(inputId = "v3cl1_1",
                         label = "Cov(V3,Cl1)", value = 0)})

      output$v3cl2_1<-shiny::renderUI({
        shiny::textInput(inputId = "v3cl2_1",
                         label = "Cov(V3,Cl2)", value = 0)})

      output$v3cl3_1<-shiny::renderUI({
        shiny::textInput(inputId = "v3cl3_1",
                         label = "Cov(V3,Cl3)", value = 0)})

      output$v1v3_1<-shiny::renderUI({
        shiny::textInput(inputId = "v1v3_1",
                         label = "Cov(V1,V3)", value = 0)})

      output$v2v3_1<-shiny::renderUI({
        shiny::textInput(inputId = "v2v3_1",
                         label = "Cov(V2,V3)", value = 0)})

      output$cl1cl3_1<-shiny::renderUI({
        shiny::textInput(inputId = "cl1cl3_1",
                         label = "Cov(Cl1,Cl3)", value = 0)})

      output$cl2cl3_1<-shiny::renderUI({
        shiny::textInput(inputId = "cl2cl3_1",
                         label = "Cov(Cl2,Cl3)", value = 0)})

      ##### Results ######
      ###### One compartment #####

      result1_1<-shiny::reactive({
        V1<-as.numeric(input$volume_mean1_1)
        Cl1<-as.numeric(input$Clearance_mean1_1)
        V1.sd<-as.numeric(input$volume_sd1_1)
        Cl1.sd<-as.numeric(input$Clearance_sd1_1)
        V1Cl1<-as.numeric(input$v1cl1_1)
        if(length(V1Cl1)==0) V1Cl1<-NA
        OneComp_Volume_Clearance(V1,Cl1,V1.sd,Cl1.sd,V1Cl1)})

      output$T1_1<-shiny::renderTable({
        Label<-c("Volume","","Clearnace","Micro Rate Constant","Exponent",
                 "Half-lives","True Coeffficient","Fractional Coefficient")
        T1.1<-data.frame(Label,result1_1())
        colnames(T1.1)[1]<-""
        T1.1},align="?",digits=4,rownames=FALSE)

      total1_1<-shiny::reactive({
        result1_1()})

      fileext <- shiny::reactive({switch(input$type1,
                                  "Excel (CSV)" = "csv",
                                  "Text (tab separated)" = "tab",
                                  "Text (Space Separated)" = "txt")})

      output$downloadData1_1<-shiny::downloadHandler(
        filename <- function() {
          paste("one_compartment_sheet1", fileext(), sep = ".")},
        content <- function(file){
          sep.t=switch(fileext(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total1_1(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE,append=FALSE)})

      ##### Two compartment #####

      result1_2<-shiny::reactive({
        V1<-as.numeric(input$volume_mean1_2)
        V2<-as.numeric(input$volume_mean2_2)
        Cl1<-as.numeric(input$Clearance_mean1_2)
        Cl2<-as.numeric(input$Clearance_mean2_2)
        V1.sd<-as.numeric(input$volume_sd1_2)
        V2.sd<-as.numeric(input$volume_sd2_2)
        Cl1.sd<-as.numeric(input$Clearance_sd1_2)
        Cl2.sd<-as.numeric(input$Clearance_sd2_2)

        V1V2<-as.numeric(input$v1v2_1)
        V1Cl1<-as.numeric(input$v1cl1_1)
        V1Cl2<-as.numeric(input$v1cl2_1)
        V2Cl1<-as.numeric(input$v2cl1_1)
        V2Cl2<-as.numeric(input$v2cl2_1)
        Cl1Cl2<-as.numeric(input$cl1cl2_1)
        TwoComp_Volume_Clearance(V1,V2,Cl1,Cl2,V1.sd,V2.sd,Cl1.sd,Cl2.sd,
                 covar=c(V1V2,V1Cl1,V1Cl2,V2Cl1,V2Cl2,Cl1Cl2))})

      output$T1_2<-shiny::renderTable({
        Label<-c("Volume","","","Clearnace","","Micro Rate Constant","","",
                 "Exponent","","Half-lives","","True Coeffficient","",
                 "Fractional Coefficient","")
        T1.2<-data.frame(Label,result1_2())
        colnames(T1.2)[1]<-""
        T1.2},align="?",digits=4,rownames=FALSE)

      total1_2<-shiny::reactive({
        result1_2()})

      output$downloadData1_2<-shiny::downloadHandler(
        filename = function() {
          paste("two_compartment_sheet1", fileext(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total1_2(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Three compartment######################

      result1_3<-shiny::reactive({
        V1<-as.numeric(input$volume_mean1_3)
        V2<-as.numeric(input$volume_mean2_3)
        V3<-as.numeric(input$volume_mean3_3)
        Cl1<-as.numeric(input$Clearance_mean1_3)
        Cl2<-as.numeric(input$Clearance_mean2_3)
        Cl3<-as.numeric(input$Clearance_mean3_3)
        V1.sd<-as.numeric(input$volume_sd1_3)
        V2.sd<-as.numeric(input$volume_sd2_3)
        V3.sd<-as.numeric(input$volume_sd3_3)
        Cl1.sd<-as.numeric(input$Clearance_sd1_3)
        Cl2.sd<-as.numeric(input$Clearance_sd2_3)
        Cl3.sd<-as.numeric(input$Clearance_sd3_3)
        V1V2<-as.numeric(input$v1v2_1);    V1V3<-as.numeric(input$v1v3_1)
        V1Cl1<-as.numeric(input$v1cl1_1);  V1Cl2<-as.numeric(input$v1cl2_1)
        V1Cl3<-as.numeric(input$v1cl3_1);  V2V3<-as.numeric(input$v2v3_1)
        V2Cl1<-as.numeric(input$v2cl1_1);  V2Cl2<-as.numeric(input$v2cl2_1)
        V2Cl3<-as.numeric(input$v2cl3_1);  V3Cl1<-as.numeric(input$v3cl1_1)
        V3Cl2<-as.numeric(input$v3cl2_1);  V3Cl3<-as.numeric(input$v3cl3_1)
        Cl1Cl2<-as.numeric(input$cl1cl2_1);Cl1Cl3<-as.numeric(input$cl1cl3_1)
        Cl2Cl3<-as.numeric(input$cl2cl3_1)
        ThreeComp_Volume_Clearance(V1,V2,V3,Cl1,Cl2,Cl3,
                 V1.sd,V2.sd,V3.sd,Cl1.sd,Cl2.sd,Cl3.sd,
                 covar=c(V1V2,V1V3,V1Cl1,V1Cl2,V1Cl3,V2V3,V2Cl1,V2Cl2,V2Cl3,
                         V3Cl1,V3Cl2,V3Cl3,Cl1Cl2,Cl1Cl3,Cl2Cl3))})

      output$T1_3<-shiny::renderTable({
        Label<-c("Volume","","","","Clearnace","","",
                 "Micro Rate Constant","","","","",
                 "Exponent","","","Half-lives","","",
                 "True Coeffficient","","","Fractional Coefficient","","")
        T1.3<-data.frame(Label,result1_3())
        colnames(T1.3)[1]<-""
        T1.3},align="?",digits=4,rownames=FALSE)

      total1_3<-shiny::reactive({
        result1_3()})

      output$downloadData1_3<-shiny::downloadHandler(
        filename = function() {
          paste("three_compartment_sheet1", fileext(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total1_3(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      #########################sheet2#########################################

      output$comp1_2_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean21_1",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd21_1",
            label = "V1 Std.err")))})

      output$comp1_2_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k10_mean_1",
            label = "k10  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k10_sd_1",
            label = "k10 Std.err")))})

      output$comp2_2_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean21_2",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd21_2",
            label = "V1 Std.err")))})

      output$comp2_2_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k10_mean_2",
            label = "k10  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k10_sd_2",
            label = "k10 Std.err")))})

      output$comp2_2_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k12_mean_2",
            label = "k12  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k12_sd_2",
            label = "k12 Std.err")))})

      output$comp2_2_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k21_mean_2",
            label = "k21  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k21_sd_2",
            label = "k21 Std.err")))})

      output$comp3_2_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean21_3",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd21_3",
            label = "V1 Std.err")))})

      output$comp3_2_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k10_mean_3",
            label = "k10  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k10_sd_3",
            label = "k10 Std.err")))})

      output$comp3_2_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k12_mean_3",
            label = "k12  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k12_sd_3",
            label = "k12 Std.err")))})

      output$comp3_2_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k13_mean_3",
            label = "k13  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k13_sd_3",
            label = "k13 Std.err")))})

      output$comp3_2_5<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k21_mean_3",
            label = "k21  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k21_sd_3",
            label = "k21 Std.err")))})

      output$comp3_2_6<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k31_mean_3",
            label = "k31  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k31_sd_3",
            label = "k31 Std.err")))})

      ###################### covariance ######################

      ######################One compartment###################

      output$v1k10_2<-shiny::renderUI({
        shiny::textInput(inputId = "v1k10_2",
                         label = "Cov(V1,k10)", value = 0)})

      ######################Two compartment######################

      output$v1k12_2<-shiny::renderUI({
        shiny::textInput(inputId = "v1k12_2",
                         label = "Cov(V1,k12)", value = 0)})

      output$v1k21_2<-shiny::renderUI({
        shiny::textInput(inputId = "v1k21_2",
                         label = "Cov(V1,k21)", value = 0)})

      output$k10k12_2<-shiny::renderUI({
        shiny::textInput(inputId = "k10k12_2",
                         label = "Cov(k10,k12)", value = 0)})

      output$k10k21_2<-shiny::renderUI({
        shiny::textInput(inputId = "k10k21_2",
                         label = "Cov(k10,k21)", value = 0)})

      output$k12k21_2<-shiny::renderUI({
        shiny::textInput(inputId = "k12k21_2",
                         label = "Cov(k12,k21)", value = 0)})

      ######################Three compartment######################

      output$v1k13_2<-shiny::renderUI({
        shiny::textInput(inputId = "v1k13_2",
                         label = "Cov(V1,k13)", value = 0)})

      output$v1k31_2<-shiny::renderUI({
        shiny::textInput(inputId = "v1k31_2",
                         label = "Cov(V1,k31)", value = 0)})

      output$k13k31_2<-shiny::renderUI({
        shiny::textInput(inputId = "k13k31_2",
                         label = "Cov(k13,k31)", value = 0)})

      output$k10k13_2<-shiny::renderUI({
        shiny::textInput(inputId = "k10k13_2",
                         label = "Cov(k10,k13)", value = 0)})

      output$k10k31_2<-shiny::renderUI({
        shiny::textInput(inputId = "k10k31_2",
                         label = "Cov(k10,k31)", value = 0)})

      output$k12k31_2<-shiny::renderUI({
        shiny::textInput(inputId = "k12k31_2",
                         label = "Cov(k12,k31)", value = 0)})

      output$k21k13_2<-shiny::renderUI({
        shiny::textInput(inputId = "k21k13_2",
                         label = "Cov(k21,k13)", value = 0)})

      output$k21k31_2<-shiny::renderUI({
        shiny::textInput(inputId = "k21k31_2",
                         label = "Cov(k21,k31)", value = 0)})

      output$k12k13_2<-shiny::renderUI({
        shiny::textInput(inputId = "k12k13_2",
                         label = "Cov(k12,k13)", value = 0)})

      ###########################################################
      #########################Results###########################
      ###########################################################

      ######################One compartment######################

      result2_1<-shiny::reactive({
        V1<-as.numeric(input$volume_mean21_1)
        k10<-as.numeric(input$k10_mean_1)
        V1.sd<-as.numeric(input$volume_sd21_1)
        k10.sd<-as.numeric(input$k10_sd_1)
        V1k10<-as.numeric(input$v1k10_2)
        OneComp_Volume_RateConstant(V1,k10,V1.sd,k10.sd,covar=c(V1k10))})

      output$T2_1<-shiny::renderTable({
        Label<-c("Volume","","Clearnace","Micro Rate Constant",
                 "Exponent","Half-lives",
                 "True Coeffficient","Fractional Coefficient")

        T2.1<-data.frame(Label,result2_1())
        colnames(T2.1)[1]<-""
        T2.1},align="?",digits=4,rownames=FALSE)

      total2_1<-shiny::reactive({result2_1()})

      fileext2 <- shiny::reactive({switch(input$type2,
                                   "Excel (CSV)" = "csv",
                                   "Text (tab separated)" = "tab",
                                   "Text (Space Separated)" = "txt")})

      output$downloadData2_1<-shiny::downloadHandler(
        filename = function() {
          paste("one_compartment_sheet2", fileext2(), sep = ".")},
        content = function(file){sep.t=switch(fileext2(),
                                              csv=",",
                                              tab="\t",
                                              txt=" ")
        utils::write.table(total2_1(), file, sep = sep.t,
                    row.names = FALSE,col.names = TRUE)})

      ######################Two compartment######################

      result2_2<-shiny::reactive({
        V1<-as.numeric(input$volume_mean21_2)
        k10<-as.numeric(input$k10_mean_2)
        k12<-as.numeric(input$k12_mean_2)
        k21<-as.numeric(input$k21_mean_2)
        V1.sd<-as.numeric(input$volume_sd21_2)
        k10.sd<-as.numeric(input$k12_sd_2)
        k12.sd<-as.numeric(input$k12_sd_2)
        k21.sd<-as.numeric(input$k21_sd_2)
        V1k10<-as.numeric(input$v1k10_2)
        V1k12<-as.numeric(input$v1k12_2)
        V1k21<-as.numeric(input$v1k21_2)
        k10k12<-as.numeric(input$k10k12_2)
        k10k21<-as.numeric(input$k10k21_2)
        k12k21<-as.numeric(input$k12k21_2)
        TwoComp_Volume_RateConstant(V1,k10,k12,k21,V1.sd,k10.sd,k12.sd,k21.sd,
                 covar=c(V1k10,V1k12,V1k21,k10k12,k10k21,k12k21))})

      output$T2_2<-shiny::renderTable({
        Label<-c("Volume","","","Clearnace","","Micro Rate Constant","","",
                 "Exponent","","Half-lives","",
                 "True Coeffficient","","Fractional Coefficient","")
        T2.2<-data.frame(Label,result2_2())
        colnames(T2.2)[1]<-""
        T2.2},align="?",digits=4,rownames=FALSE)

      total2_2<-shiny::reactive({result2_2()})

      output$downloadData2_2<-shiny::downloadHandler(
        filename = function() {
          paste("two_compartment_sheet2", fileext2(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext2(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total2_2(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Three compartment######################

      result2_3<-shiny::reactive({
        V1<-as.numeric(input$volume_mean21_3)
        k10<-as.numeric(input$k10_mean_3)
        k12<-as.numeric(input$k12_mean_3)
        k13<-as.numeric(input$k13_mean_3)
        k21<-as.numeric(input$k21_mean_3)
        k31<-as.numeric(input$k31_mean_3)
        V1.sd<-as.numeric(input$volume_sd21_3)
        k10.sd<-as.numeric(input$k10_sd_3);
        k12.sd<-as.numeric(input$k12_sd_3)
        k13.sd<-as.numeric(input$k13_sd_3);
        k21.sd<-as.numeric(input$k21_sd_3)
        k31.sd<-as.numeric(input$k31_sd_3)

        V1k10<-as.numeric(input$v1k10_2);
        V1k12<-as.numeric(input$v1k12_2)
        V1k13<-as.numeric(input$v1k13_2);
        V1k21<-as.numeric(input$v1k21_2)
        V1k31<-as.numeric(input$v1k31_2);
        k10k12<-as.numeric(input$k10k12_2)
        k10k13<-as.numeric(input$k10k13_2);
        k10k21<-as.numeric(input$k10k21_2)
        k10k31<-as.numeric(input$k10k31_2);
        k12k13<-as.numeric(input$k12k13_2)
        k12k21<-as.numeric(input$k12k21_2);
        k12k31<-as.numeric(input$k12k31_2)
        k21k31<-as.numeric(input$k21k31_2);
        k13k21<-as.numeric(input$k13k21_2)
        k13k31<-as.numeric(input$k13k31_2)

        ThreeComp_Volume_RateConstant(V1,k10,k12,k13,k21,k31,V1.sd,
          k10.sd,k12.sd,k13.sd,k21.sd,k31.sd,
          covar=c(V1k10,V1k12,V1k13,V1k21,V1k31,k10k12,k10k13,k10k21,k10k31,
            k12k13,k12k21,k12k31,k13k21,k13k31,k21k31))})

      output$T2_3<-shiny::renderTable({
        Label<-c("Volume","","","","Clearnace","","",
                 "Micro Rate Constant","","","","",
                 "Exponent","","","Half-lives","","",
                 "True Coeffficient","","","Fractional Coefficient","","")
        T2.3<-data.frame(Label,result2_3())
        colnames(T2.3)[1]<-""
        T2.3},align="?",digits=4,rownames=FALSE)

      total2_3<-shiny::reactive({result2_3()})

      output$downloadData2_3<-shiny::downloadHandler(
        filename = function() {
          paste("three_compartment_sheet2", fileext2(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext2(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total2_3(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################sheet3############################################

      output$comp1_3_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Cl1_mean_1",
            label = "Cl1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Cl1_sd_1",
            label = "Cl1 Std.err")))})

      output$comp1_3_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "t_alpha_mean_1",
            label = "t_alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "t_alpha_sd_1",
            label = "t_alpha Std.err")))})

      output$comp2_3_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean31_2",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd31_2",
            label = "V1 Std.err")))})

      output$comp2_3_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Cl1_mean_2",
            label = "Cl1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Cl1_sd_2",
            label = "Cl1 Std.err")))})

      output$comp2_3_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "t_alpha_mean_2",
            label = "t_alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "t_alpha_sd_2",
            label = "t_alpha Std.err")))})

      output$comp2_3_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "t_beta_mean_2",
            label = "t_beta  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "t_beta_sd_2",
            label = "t_beta Std.err")))})

      output$comp3_3_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean31_3",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd31_3",
            label = "V1 Std.err")))})

      output$comp3_3_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Vd_mean_3",
            label = "Vd  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Vd_sd_3",
            label = "Vd Std.err")))})

      output$comp3_3_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "Cl1_mean_3",
            label = "Cl1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "Cl1_sd_3",
            label = "Cl1 Std.err")))})

      output$comp3_3_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "t_alpha_mean_3",
            label = "t_alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "t_alpha_sd_3",
            label = "t_alpha Std.err")))})

      output$comp3_3_5<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "t_beta_mean_3",
            label = "t_beta  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "t_beta_sd_3",
            label = "t_beta Std.err")))})

      output$comp3_3_6<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "t_gamma_mean_3",
            label = "t_gamma  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "t_gamma_sd_3",
            label = "t_gamma Std.err")))})


      ###################### covariance ######################
      ######################One compartment######################

      output$cl1talpha_3<-shiny::renderUI({
        shiny::textInput(inputId = "cl1talpha_3",
                  label = "Cov(Cl1,t_alpha)", value = 0) })

      ######################Two compartment######################

      output$v1cl1_3<-shiny::renderUI({
        shiny::textInput(inputId = "v1cl1_3",
                  label = "Cov(V1,Cl1)", value = 0)})

      output$cl1tbeta_3<-shiny::renderUI({
        shiny::textInput(inputId = "cl1tbeta_3",
                  label = "Cov(Cl1,t_beta)", value = 0)})

      output$v1talpha_3<-shiny::renderUI({
        shiny::textInput(inputId = "v1talpha_3",
                  label = "Cov(V1,t_alpha)", value = 0)})

      output$v1tbeta_3<-shiny::renderUI({
        shiny::textInput(inputId = "v1tbeta_3",
                  label = "Cov(V1,t_beta)", value = 0)})

      output$talphatbeta_3<-shiny::renderUI({
        shiny::textInput(inputId = "talphatbeta_3",
                  label = "Cov(t_alpha,t_beta)",value = 0)})

      ######################Three compartment######################

      output$v1vd_3<-shiny::renderUI({
        shiny::textInput(inputId = "v1vd_3",
                  label = "Cov(V1,Vd)", value = 0)})

      output$v1tgamma_3<-shiny::renderUI({
        shiny::textInput(inputId = "v1tgamma_3",
                  label = "Cov(V1,t_gamma)", value = 0)})

      output$vdtgamma_3<-shiny::renderUI({
        shiny::textInput(inputId = "vdtgamma_3",
                  label = "Cov(Vd,t_gamma)", value = 0)})

      output$cl1tgamma_3<-shiny::renderUI({
        shiny::textInput(inputId = "cl1tgamma_3",
                  label = "Cov(Cl1,t_gamma)",value = 0)})

      output$talphatgamma_3<-shiny::renderUI({
        shiny::textInput(inputId = "talphatgamma_3",
                  label = "Cov(t_alpha,t_gamma)", value = 0)})

      output$tbetatgamma_3<-shiny::renderUI({
        shiny::textInput(inputId = "tbetatgamma_3",
                  label = "Cov(t_beta,t_gamma)",value = 0)})

      output$vdtalpha_3<-shiny::renderUI({
        shiny::textInput(inputId = "vdtalpha_3",
                  label = "Cov(Vd,t_alpha)", value = 0)})

      output$vdtbeta_3<-shiny::renderUI({
        shiny::textInput(inputId = "vdtbeta_3",
                  label = "Cov(Vd,t_beta)", value = 0)})

      output$vdcl1_3<-shiny::renderUI({
        shiny::textInput(inputId = "vdcl1_3",
                         label = "Cov(Vd,Cl1)", value = 0)})

      ###########################################################
      #########################Results###########################
      ###########################################################

      ######################One compartment######################

      result3_1<-shiny::reactive({
        Cl1<-as.numeric(input$Cl1_mean_1)
        t_alpha<-as.numeric(input$t_alpha_mean_1)
        Cl1.sd<-as.numeric(input$Cl1_sd_1)
        t_alpha.sd<-as.numeric(input$t_alpha_sd_1)
        Cl1talpha<-as.numeric(input$cl1talpha_3)
        OneComp_Volume_Clearance_HalfLife(Cl1,t_alpha,Cl1.sd,t_alpha.sd,
           covar=c(Cl1talpha))})

      output$T3_1<-shiny::renderTable({
        Label<-c("Volume","","Clearnace","Micro Rate Constant",
                 "Exponent","Half-lives",
                 "True Coeffficient","Fractional Coefficient")

        T3.1<-data.frame(Label,result3_1())
        colnames(T3.1)[1]<-""
        T3.1},align="?",digits=4,rownames=FALSE)

      total3_1<-shiny::reactive({result3_1() })

      fileext3 <- shiny::reactive({switch(input$type3,
                                   "Excel (CSV)" = "csv",
                                   "Text (tab separated)" = "tab",
                                   "Text (Space Separated)" = "txt")})

      output$downloadData3_1<-shiny::downloadHandler(
        filename = function() {
          paste("one_compartment_sheet3", fileext3(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext3(),
                       csv=",",
                       tab="\t",
                       txt=" ")
         utils::write.table(total3_1(), file, sep = sep.t,
                        row.names = FALSE,col.names = TRUE)})

      ####################Two compartment######################

      result3_2<-shiny::reactive({
        V1<-as.numeric(input$volume_mean31_2)
        Cl1<-as.numeric(input$Cl1_mean_2)
        t_alpha<-as.numeric(input$t_alpha_mean_2)
        t_beta<-as.numeric(input$t_beta_mean_2)
        V1.sd<-as.numeric(input$volume_sd31_2)
        Cl1.sd<-as.numeric(input$Cl1_sd_2)
        t_alpha.sd<-as.numeric(input$t_alpha_sd_2)
        t_beta.sd<-as.numeric(input$t_beta_sd_2)

        V1Cl1<-as.numeric(input$v1cl1_2)
        V1talpha<-as.numeric(input$v1talpha_2)
        V1tbeta<-as.numeric(input$v1tbeta_2)
        Cl1talpha<-as.numeric(input$cl1talpha_2)
        Cl1tbeta<-as.numeric(input$cl1tbeta_2)
        talphatbeta<-as.numeric(input$talphatbeta_2)

        TwoComp_Volume_Clearance_HalfLife(V1,Cl1,t_alpha,t_beta,
          V1.sd,Cl1.sd,t_alpha.sd,t_beta.sd,
          covar=c(V1Cl1,V1talpha,V1tbeta,Cl1talpha,Cl1tbeta,talphatbeta))})

      output$T3_2<-shiny::renderTable({
        Label<-c("Volume","","","Clearnace","",
                 "Micro Rate Constant","","",
                 "Exponent","","Half-lives","",
                 "True Coeffficient","","Fractional Coefficient","")
        T3.2<-data.frame(Label,result3_2())
        colnames(T3.2)[1]<-""
        T3.2},align="?",digits=4,rownames=FALSE)

      total3_2<-shiny::reactive({ result3_2()})

      output$downloadData3_2<-shiny::downloadHandler(
        filename = function() {
          paste("two_compartment_sheet3", fileext3(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext3(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total3_2(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Three compartment######################

      result3_3<-shiny::reactive({
        V1<-as.numeric(input$volume_mean31_3)
        Vd<-as.numeric(input$Vd_mean_3)
        Cl1<-as.numeric(input$Cl1_mean_3)
        t_alpha<-as.numeric(input$t_alpha_mean_3)
        t_beta<-as.numeric(input$t_beta_mean_3)
        t_gamma<-as.numeric(input$t_gamma_mean_3)

        V1.sd<-as.numeric(input$volume_sd31_3)
        Vd.sd<-as.numeric(input$Vd_sd_3)
        Cl1.sd<-as.numeric(input$Cl1_sd_3)
        t_alpha.sd<-as.numeric(input$t_alpha_sd_3)
        t_beta.sd<-as.numeric(input$t_beta_sd_3)
        t_gamma.sd<-as.numeric(input$t_gamma_sd_3)

        V1Vd<-as.numeric(input$v1vd_3)
        V1Cl1<-as.numeric(input$v1cl1_3)
        V1talpha<-as.numeric(input$v1talpha_3)
        V1tbeta<-as.numeric(input$v1tbeta_3)
        V1tgamma<-as.numeric(input$v1tgamma_3)
        VdCl1<-as.numeric(input$vdcl1_3)
        Vdtalpha<-as.numeric(input$vdtalpha_3)
        Vdtbeta<-as.numeric(input$vdtbeta_3)
        Vdtgamma<-as.numeric(input$vdtgamma_3)
        Cl1talpha<-as.numeric(input$cl1talpha_3)
        Cl1tbeta<-as.numeric(input$cl1tbeta_3)
        Cl1tgamma<-as.numeric(input$cl1tgamma_3)
        talphatbeta<-as.numeric(input$talphatbeta_3)
        talphatgamma<-as.numeric(input$talphatgamma_3)
        tbetatgamma<-as.numeric(input$tbetatgamma_3)
        ThreeComp_Volume_Clearance_HalfLife(V1,Vd,Cl1,t_alpha,t_beta,t_gamma,
          V1.sd,Vd.sd,Cl1.sd,t_alpha.sd,t_beta.sd,t_gamma.sd,
          covar=c(V1Vd,V1Cl1,V1talpha,V1tbeta,V1tgamma,
            VdCl1,Vdtalpha,Vdtbeta,Vdtgamma,Cl1talpha,Cl1tbeta,Cl1tgamma,
            talphatbeta,talphatgamma,tbetatgamma))})

      output$T3_3<-shiny::renderTable({
        Label<-c("Volume","","","","Clearnace","","",
                 "Micro Rate Constant","","","","",
                 "Exponent","","","Half-lives","","",
                 "True Coeffficient","","","Fractional Coefficient","","")
        T3.3<-data.frame(Label,result3_3())
        colnames(T3.3)[1]<-""
        T3.3},align="?",digits=4,rownames=FALSE)

      total3_3<-shiny::reactive({result3_3()})

      output$downloadData3_3<-shiny::downloadHandler(
        filename = function() {
          paste("three_compartment_sheet3", fileext3(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext3(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total3_3(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ##############sheet4############################################

      output$comp1_4_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "A_mean_1",
            label = "A  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "A_sd_1",
            label = "A Std.err")))})

      output$comp1_4_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "alpha_mean_1",
            label = "alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "alpha_sd_1",
            label = "alpha Std.err")))})

      output$comp1_4_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "A_mean_1",
            label = "A  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "A_sd_1",
            label = "A Std.err")))})

      output$comp2_4_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "A_mean_2",
            label = "A  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "A_sd_2",
            label = "A Std.err")))})

      output$comp2_4_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "B_mean_2",
            label = "B  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "B_sd_2",
            label = "B Std.err")))})

      output$comp2_4_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "A_mean_2",
            label = "A  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "A_sd_2",
            label = "A Std.err")))})

      output$comp2_4_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "B_mean_2",
            label = "B  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "B_sd_2",
            label = "B Std.err")))})

      output$comp2_4_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "alpha_mean_2",
            label = "alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "alpha_sd_2",
            label = "alpha Std.err")))})

      output$comp2_4_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "beta_mean_2",
            label = "beta  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "beta_sd_2",
            label = "beta Std.err")))})

      output$comp3_4_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "A_mean_3",
            label = "A  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "A_sd_3",
            label = "A Std.err")))})

      output$comp3_4_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "B_mean_3",
            label = "B  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "B_sd_3",
            label = "B Std.err")))})

      output$comp3_4_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "C_mean_3",
            label = "C  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "C_sd_3",
            label = "C Std.err")))})

      output$comp3_4_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "alpha_mean_3",
            label = "alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "alpha_sd_3",
            label = "alpha Std.err")))})

      output$comp3_4_5<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "beta_mean_3",
            label = "beta  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "beta_sd_3",
            label = "beta Std.err")))})

      output$comp3_4_6<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "gamma_mean_3",
            label = "gamma  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "gamma_sd_3",
            label = "gamma Std.err")))})

      ###################### covariance ######################

      ######################One compartment######################

      output$aalpha_4<-shiny::renderUI({
        shiny::textInput(inputId = "aalpha_4",
                  label = "Cov(A,alpha)", value = 0)})

      ######################Two compartment######################

      output$ab_4<-shiny::renderUI({
        shiny::textInput(inputId = "ab_4",
                  label = "Cov(A,B)", value = 0)})

      output$abeta_4<-shiny::renderUI({
        shiny::textInput(inputId = "abeta_4",
                  label = "Cov(A,beta)", value = 0)})

      output$balpha_4<-shiny::renderUI({
        shiny::textInput(inputId = "balpha_4",
                  label = "Cov(B,alpha)", value = 0)})

      output$bbeta_4<-shiny::renderUI({
        shiny::textInput(inputId = "bbeta_4",
                  label = "Cov(B,beta)", value = 0)})

      output$alphabeta_4<-shiny::renderUI({
        shiny::textInput(inputId = "alphabeta_4",
                  label = "Cov(alpha,beta)", value = 0)})

      ######################Three compartment######################

      output$ac_4<-shiny::renderUI({
        shiny::textInput(inputId = "ac_4",
                  label = "Cov(A,C)", value = 0)})

      output$agamma_4<-shiny::renderUI({
        shiny::textInput(inputId = "agamma_4",
                  label = "Cov(A,gamma)", value = 0)})

      output$bgamma_4<-shiny::renderUI({
        shiny::textInput(inputId = "bgamma_4",
                  label = "Cov(B,gamma)", value = 0)})

      output$cbeta_4<-shiny::renderUI({
        shiny::textInput(inputId = "cbeta_4",
                  label = "Cov(C,beta)", value = 0)})

      output$cgamma_4<-shiny::renderUI({
        shiny::textInput(inputId = "cgamma_4",
                  label = "Cov(C,gamma)", value = 0)})

      output$alphagamma_4<-shiny::renderUI({
        shiny::textInput(inputId = "alphagamma_4",
                  label = "Cov(alpha,gamma)", value = 0)})

      output$betagamma_4<-shiny::renderUI({
        shiny::textInput(inputId = "betagamma_4",
                  label = "Cov(beta,gamma)", value = 0)})

      output$bc_4<-shiny::renderUI({
        shiny::textInput(inputId = "bc_4",
                  label = "Cov(B,C)", value = 0)})

      output$calpha_4<-shiny::renderUI({
        shiny::textInput(inputId = "calpha_4",
                  label = "Cov(C,alpha)", value = 0)})

      ###########################################################
      #########################Results###########################
      ###########################################################

      ######################One compartment######################

      result4_1<-shiny::reactive({
        A<-as.numeric(input$A_mean_1)
        alpha<-as.numeric(input$alpha_mean_1)
        A.sd<-as.numeric(input$A_sd_1)
        alpha.sd<-as.numeric(input$alpha_sd_1)
        Aalpha<-as.numeric(input$aalpha_4)
        OneComp_Coefficient_Exponent(A,alpha,A.sd,alpha.sd,covar=Aalpha)})

      output$T4_1<-shiny::renderTable({
        Label<-c("Volume","","Clearnace","Micro Rate Constant",
                 "Exponent","Half-lives",
               "True Coeffficient","Fractional Coefficient")

        T4.1<-data.frame(Label,result4_1())
        colnames(T4.1)[1]<-""
        T4.1},align="?",digits=4,rownames=FALSE)


      total4_1<-shiny::reactive({result4_1()})

      fileext4 <- shiny::reactive({switch(input$type4,
                                   "Excel (CSV)" = "csv",
                                   "Text (tab separated)" = "tab",
                                   "Text (Space Separated)" = "txt")})

      output$downloadData4_1<-shiny::downloadHandler(
        filename = function() {
          paste("one_compartment_sheet4", fileext4(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext4(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total4_1(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Two compartment######################

      result4_2<-shiny::reactive({
        A<-as.numeric(input$A_mean_2)
        B<-as.numeric(input$B_mean_2)
        alpha<-as.numeric(input$alpha_mean_2)
        beta<-as.numeric(input$beta_mean_2)
        A.sd<-as.numeric(input$A_sd_2)
        B.sd<-as.numeric(input$B_sd_2)
        alpha.sd<-as.numeric(input$alpha_sd_2)
        beta.sd<-as.numeric(input$beta_sd_2)
        AB<-as.numeric(input$ab_4)
        Aalpha<-as.numeric(input$aalpha_4)
        Abeta<-as.numeric(input$abeta_4)
        Balpha<-as.numeric(input$balpha_4)
        Bbeta<-as.numeric(input$bbeta_4)
        alphabeta<-as.numeric(input$alphabeta_4)
        TwoComp_Coefficient_Exponent(A,B,alpha,beta,A.sd,B.sd,alpha.sd,beta.sd,
                 covar=c(AB,Aalpha,Abeta,Balpha,Bbeta,alphabeta))})

      output$T4_2<-shiny::renderTable({
        Label<-c("Volume","","","Clearnace","",
                 "Micro Rate Constant","","",
                 "Exponent","","Half-lives","",
                 "True Coeffficient","","Fractional Coefficient","")
        T4.2<-data.frame(Label,result4_2())
        colnames(T4.2)[1]<-""
        T4.2},align="?",digits=4,rownames=FALSE)

      total4_2<-shiny::reactive({result4_2()})

      output$downloadData4_2<-shiny::downloadHandler(
        filename = function() {
          paste("two_compartment_sheet4", fileext4(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext4(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total4_2(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Three compartment######################

      result4_3<-shiny::reactive({
        A<-as.numeric(input$A_mean_3)
        B<-as.numeric(input$B_mean_3)
        C<-as.numeric(input$C_mean_3)
        alpha<-as.numeric(input$alpha_mean_3)
        beta<-as.numeric(input$beta_mean_3)
        gamma<-as.numeric(input$gamma_mean_3)
        A.sd<-as.numeric(input$A_sd_3)
        B.sd<-as.numeric(input$B_sd_3)
        C.sd<-as.numeric(input$C_sd_3)
        alpha.sd<-as.numeric(input$alpha_sd_3)
        beta.sd<-as.numeric(input$beta_sd_3)
        gamma.sd<-as.numeric(input$gamma_sd_3)

        AB<-as.numeric(input$ab_4);
        AC<-as.numeric(input$ac_4)
        Aalpha<-as.numeric(input$aalpha_4);
        Abeta<-as.numeric(input$abeta_4)
        Agamma<-as.numeric(input$agamma_4);
        BC<-as.numeric(input$bc_4)
        Balpha<-as.numeric(input$balpha_4);
        Bbeta<-as.numeric(input$bbeta_4)
        Bgamma<-as.numeric(input$bgamma_4);
        Calpha<-as.numeric(input$calpha_4)
        Cbeta<-as.numeric(input$cbeta_4);
        Cgamma<-as.numeric(input$cgamma_4)
        alphabeta<-as.numeric(input$alphabeta_4);
        alphagamma<-as.numeric(input$alphagamma_4)
        betagamma<-as.numeric(input$betagamma_4);

        ThreeComp_Coefficient_Exponent(A,B,C,alpha,beta,gamma,A.sd,B.sd,C.sd,
          alpha.sd,beta.sd,gamma.sd,
          covar=c(AB,AC,Aalpha,Abeta,Agamma,BC,Balpha,Bbeta,Bgamma,Calpha,
            Cbeta, Cgamma,alphabeta,alphagamma,betagamma))})

      output$T4_3<-shiny::renderTable({
        Label<-c("Volume","","","","Clearnace","","",
                 "Micro Rate Constant","","","","",
                 "Exponent","","","Half-lives","","",
                 "True Coeffficient","","","Fractional Coefficient","","")
        T4.3<-data.frame(Label,result4_3())
        colnames(T4.3)[1]<-""
        T4.3},align="?",digits=4,rownames=FALSE)

      total4_3<-shiny::reactive({result4_3()})

      output$downloadData4_3<-shiny::downloadHandler(
        filename = function() {
          paste("three_compartment_sheet4", fileext4(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext4(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total4_3(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      #####################sheet5############################################

      output$comp1_5_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean51_1",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd51_1",
            label = "V1 Std.err")))})

      output$comp1_5_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "alpha_mean51_1",
            label = "alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "alpha_sd51_1",
            label = "alpha Std.err")))})

      output$comp2_5_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean51_2",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd51_2",
            label = "V1 Std.err")))})

      output$comp2_5_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "alpha_mean51_2",
            label = "alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "alpha_sd51_2",
            label = "alpha Std.err")))})

      output$comp2_5_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "beta_mean51_2",
            label = "beta  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "beta_sd51_2",
            label = "beta Std.err")))})

      output$comp2_5_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k21_mean51_2",
            label = "k21  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k21_sd51_2",
            label = "k21 Std.err")))})

      output$comp3_5_1<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "volume_mean51_3",
            label = "V1  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "volume_sd51_3",
            label = "V1 Std.err")))})

      output$comp3_5_2<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "alpha_mean51_3",
            label = "alpha  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "alpha_sd51_3",
            label = "alpha Std.err")))})

      output$comp3_5_3<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "beta_mean51_3",
            label = "beta  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "beta_sd51_3",
            label = "beta Std.err")))})

      output$comp3_5_4<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "gamma_mean51_3",
            label = "gamma  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "gamma_sd51_3",
            label = "gamma Std.err")))})

      output$comp3_5_5<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k21_mean51_3",
            label = "k21  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k21_sd51_3",
            label = "k21 Std.err")))})

      output$comp3_5_6<-shiny::renderUI({
        shiny::fluidRow(shiny::column(6,
          shiny::textInput(inputId = "k31_mean51_3",
            label = "k31  Estimate")),
          shiny::column(6,shiny::textInput(inputId = "k31_sd51_3",
            label = "k31 Std.err")))})


      ###################### covariance ######################
      ######################One compartment######################

      output$v1alpha_5<-shiny::renderUI({
        shiny::textInput(inputId = "v1alpha_5",
                  label = "Cov(V1,alpha)", value = 0)})

      ######################Two compartment######################

      output$v1beta_5<-shiny::renderUI({
        shiny::textInput(inputId = "v1beta_5",
                  label = "Cov(V1,beta)", value = 0)})

      output$v1k21_5<-shiny::renderUI({
        shiny::textInput(inputId = "v1k21_5",
                  label = "Cov(V1,k21)", value = 0)})

      output$alphabeta_5<-shiny::renderUI({
        shiny::textInput(inputId = "alphabeta_5",
                  label = "Cov(alpha,beta)", value = 0)})

      output$alphak21_5<-shiny::renderUI({
        shiny::textInput(inputId = "alphak21_5",
                  label = "Cov(alpha,k21)", value = 0)})

      output$betak21_5<-shiny::renderUI({
        shiny::textInput(inputId = "betak21_5",
                  label = "Cov(beta,k21)", value = 0)})

      ######################Three compartment######################

      output$v1gamma_5<-shiny::renderUI({
        shiny::textInput(inputId = "v1gamma_5",
                  label = "Cov(V1,gamma)", value = 0)})

      output$v1k31_5<-shiny::renderUI({
        shiny::textInput(inputId = "v1k31_5",
                  label = "Cov(V1,k31)", value = 0)})

      output$alphagamma_5<-shiny::renderUI({
        shiny::textInput(inputId = "alphagamma_5",
                  label = "Cov(alpha,gamma)", value = 0)})

      output$alphak31_5<-shiny::renderUI({
        shiny::textInput(inputId = "alphak31_5",
                  label = "Cov(alpha,k31)", value = 0)})

      output$betagamma_5<-shiny::renderUI({
        shiny::textInput(inputId = "betagamma_5",
                  label = "Cov(beta,gamma)", value = 0)})

      output$betak31_5<-shiny::renderUI({
        shiny::textInput(inputId = "betak31_5",
                  label = "Cov(beta,k31)", value = 0)})

      output$gammak21_5<-shiny::renderUI({
        shiny::textInput(inputId = "gammak21_5",
                  label = "Cov(gamma,k21)", value = 0)})

      output$gammak31_5<-shiny::renderUI({
        shiny::textInput(inputId = "gammak31_5",
                  label="Cov(gamma,k31)", value=0)})

      output$k21k31_5<-shiny::renderUI({
        shiny::textInput(inputId = "k21k31_5",
                  label="Cov(k21,k31)", value = 0)})


      ###########################################################
      #########################Results###########################
      ###########################################################

      ######################One compartment######################

      result5_1<-shiny::reactive({
        V1<-as.numeric(input$volume_mean51_1)
        V1.sd<-as.numeric(input$volume_sd51_1)
        alpha<-as.numeric(input$alpha_mean51_1)
        alpha.sd<-as.numeric(input$alpha_sd51_1)
        V1alpha<-as.numeric(input$v1alpha_5)
        OneComp_Volume_Exponent(V1,alpha,V1.sd,alpha.sd,covar=c(V1alpha))})

      output$T5_1<-shiny::renderTable({
        Label<-c("Volume","","Clearnace","Micro Rate Constant",
                 "Exponent","Half-lives",
                 "True Coeffficient","Fractional Coefficient")

        T5.1<-data.frame(Label,result5_1())
        colnames(T5.1)[1]<-""
        T5.1},align="?",digits=4,rownames=FALSE)

      total5_1<-shiny::reactive({result5_1()})

      fileext5 <- shiny::reactive({switch(input$type5,
                                   "Excel (CSV)" = "csv",
                                   "Text (tab separated)" = "tab",
                                   "Text (Space Separated)" = "txt")})

      output$downloadData5_1<-shiny::downloadHandler(
        filename = function() {
          paste("one_compartment_sheet5", fileext5(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext5(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total5_1(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Two compartment######################

      result5_2<-shiny::reactive({
        V1<-as.numeric(input$volume_mean51_2)
        alpha<-as.numeric(input$alpha_mean51_2)
        beta<-as.numeric(input$beta_mean51_2)
        k21<-as.numeric(input$k21_mean51_2)
        V1.sd<-as.numeric(input$volume_sd51_2)
        alpha.sd<-as.numeric(input$alpha_sd51_2)
        beta.sd<-as.numeric(input$beta_sd51_2)
        k21.sd<-as.numeric(input$k21_sd51_2)

        V1alpha<-as.numeric(input$v1alpha_5)
        V1beta<-as.numeric(input$v1beta_5)
        V1k21<-as.numeric(input$v1k21_5)
        alphabeta<-as.numeric(input$alphabeta_5)
        alphak21<-as.numeric(input$alphak21_5)
        betak21<-as.numeric(input$betak21_5)
        TwoComp_Volume_Exponent(V1,alpha,beta,k21,V1.sd,alpha.sd,beta.sd,k21.sd,
                 covar=c(V1alpha,V1beta,V1k21,alphabeta,alphak21,betak21)) })

      output$T5_2<-shiny::renderTable({
        Label<-c("Volume","","","Clearnace","",
                 "Micro Rate Constant","","",
                 "Exponent","","Half-lives","",
                 "True Coeffficient","","Fractional Coefficient","")
        T5.2<-data.frame(Label,result5_2())
        colnames(T5.2)[1]<-""
        T5.2},align="?",digits=4,rownames=FALSE)

      total5_2<-shiny::reactive({ result5_2()})

      output$downloadData5_2<-shiny::downloadHandler(
        filename = function() {
          paste("two_compartment_sheet5", fileext5(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext5(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total5_2(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      ######################Three compartment######################

      result5_3<-shiny::reactive({
        V1<-as.numeric(input$volume_mean51_3)
        alpha<-as.numeric(input$alpha_mean51_3)
        beta<-as.numeric(input$beta_mean51_3)
        gamma<-as.numeric(input$gamma_mean51_3)
        k21<-as.numeric(input$k21_mean51_3)
        k31<-as.numeric(input$k31_mean51_3)
        V1.sd<-as.numeric(input$volume_sd51_3)
        alpha.sd<-as.numeric(input$alpha_sd51_3)
        beta.sd<-as.numeric(input$beta_sd51_3)
        gamma.sd<-as.numeric(input$gamma_sd51_3)
        k21.sd<-as.numeric(input$k21_sd51_3)
        k31.sd<-as.numeric(input$k31_sd51_3)

        V1alpha<-as.numeric(input$v1alpha_5);
        V1beta<-as.numeric(input$v1beta_5)
        V1gamma<-as.numeric(input$v1gamma_5);
        V1k21<-as.numeric(input$v1k21_5)
        V1k31<-as.numeric(input$v1k31_5);
        alphabeta<-as.numeric(input$alphabeta_5)
        alphagamma<-as.numeric(input$alphagamma_5)
        alphak21<-as.numeric(input$alphak21_5)
        alphak31<-as.numeric(input$alphak31_5)
        betagamma<-as.numeric(input$betagamma_5)
        betak21<-as.numeric(input$betak21_5)
        betak31<-as.numeric(input$betak31_5)
        gammak21<-as.numeric(input$gammak21_5)
        gammak31<-as.numeric(input$gammak31_5)
        k21k31<-as.numeric(input$k21k31_5)
        ThreeComp_Volume_Exponent(V1,alpha,beta,gamma,k21,k31,
          V1.sd,alpha.sd,beta.sd,gamma.sd,
          k21.sd,k31.sd,covar=c(V1alpha,V1beta,V1gamma,V1k21,V1k31,alphabeta,
            alphagamma,alphak21,alphak31,betagamma,betak21,betak31,gammak21,
            gammak31,k21k31))})

      output$T5_3<-shiny::renderTable({
        Label<-c("Volume","","","","Clearnace","","",
                 "Micro Rate Constant","","","","",
                 "Exponent","","","Half-lives","","",
                 "True Coeffficient","","","Fractional Coefficient","","")
        T5.3<-data.frame(Label,result3_3())
        colnames(T5.3)[1]<-""
        T5.3},align="?",digits=4,rownames=FALSE)

      total5_3<-shiny::reactive({ result5_3()})

      output$downloadData5_3<-shiny::downloadHandler(
        filename = function() {
          paste("three_compartment_sheet5", fileext5(), sep = ".")},
        content = function(file){
          sep.t=switch(fileext5(),
                       csv=",",
                       tab="\t",
                       txt=" ")
          utils::write.table(total5_3(), file, sep = sep.t,
                      row.names = FALSE,col.names = TRUE)})

      Dataset<-shiny::reactive({
        if(is.null(input$inputfile)){return(data.frame())}
        outputData<-data.frame(do.call("read.csv",
          list(input$inputfile$datapath)))
        return(outputData)})

      output$choose_PK<-shiny::renderUI({
        if(is.null(input$inputfile))
          return()
        if(identical(Dataset(),'')||identical(Dataset(),data.frame()))
          return(NULL)
        PKconverter.env$outputData<-Dataset()
        PKconverter.env$labelList<-switch(
          paste(input$modelF1,input$modelF2,sep="-"),
          "one-M1"=c("V1","Cl1"),
          "one-M2"=c("V1","k10"),
          "one-M3"=c("Cl1","t_alpha"),
          "one-M4"=c("True_A","alpha"),
          "one-M5"=c("V1","alpha"),
          "two-M1"=c("V1","V2","Cl1","Cl2"),
          "two-M2"=c("V1","k10","k12","k21"),
          "two-M3"=c("V1","Cl1","t_alpha","t_beta"),
          "two-M4"=c("True_A","True_B","alpha","beta"),
          "two-M5"=c("V1","k21","alpha","beta"),
          "three-M1"=c("V1","V2","V3","Cl1","Cl2","Cl3"),
          "three-M2"=c("V1","k10","k12","k21","k13","k31"),
          "three-M3"=c("V1","Cl1","t_alpha","t_beta","t_gamma"),
          "three-M4"=c("True_A","True_B","True_C","alpha","beta","gamma"),
          "three-M5"=c("V1","k21","k31","alpha","beta","gamma"))
        PKconverter.env$labelList<-c("ID",PKconverter.env$labelList)

        MakeSelectList<-function(X=1){
          shiny::selectInput(paste0("PK",as.character(X)),
                      label=PKconverter.env$labelList[X],
                      choices=c("",colnames(PKconverter.env$outputData)),
                      selected="")
        }
        WidgetVector<-shiny::reactive({
          n<-length(PKconverter.env$labelList)
          lapply(X=1:n,FUN=MakeSelectList)
        })

        output$choose_PK<-shiny::renderUI({
          if(is.null(input$inputfile)|is.null(PKconverter.env$outputData)){
            return()
          } else{
            shiny::tagList(WidgetVector())}})})

      output$ShowResult<-shiny::renderTable({
        if(is.null(input$inputfile)) return()
        if(identical(Dataset(),'')||identical(Dataset(),data.frame()))
          return(NULL)
        if(length(PKconverter.env$labelList)==0) return()
        if(length(input[["PK1"]])==0) return()
        n<-length(PKconverter.env$labelList)
        if(input[[paste0("PK",n)]]=="") return()
        param.list<-rep("",n)
        for(i in 1:n){
          param.list[i]<-input[[paste0("PK",i)]]
        }
        PKconverter.env$param.list<-param.list
        sel.data<-PKconverter.env$outputData[,param.list]

        colnames(sel.data)<-PKconverter.env$labelList
        modelF<-paste(input$modelF1,input$modelF2,sep="-")
        tot.result<-NULL
        if(input$modelF1=="one"){
          for(i in 1:nrow(sel.data)){
            if(input$modelF2=="M1"){
              temp.result<-OneComp_Volume_Clearance(sel.data[i,2],
                                                    sel.data[i,3])
            } else if(input$modelF2=="M2"){
              temp.result<-OneComp_Volume_RateConstant(sel.data[i,2],
                                                       sel.data[i,3])
            } else if(input$modelF2=="M3"){
              temp.result<-OneComp_Volume_Clearance_HalfLife(sel.data[i,2],
                                                             sel.data[i,3])
            } else if(input$modelF2=="M4"){
              temp.result<-OneComp_Coefficient_Exponent(sel.data[i,2],
                                                        sel.data[i,3])
            } else if(input$modelF2=="M5"){
              temp.result<-OneComp_Volume_Exponent(sel.data[i,2],
                                                   sel.data[i,3])
            }
            tot.result<-rbind(tot.result,c(i,temp.result[,2]))
            if(i==1) colnames(tot.result)<-c("ID",rownames(temp.result))
          }
        } else if(input$modelF1=="two"){
          for(i in 1:nrow(sel.data)){
            if(input$modelF2=="M1"){
              temp.result<-TwoComp_Volume_Clearance(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],sel.data[i,5])
            } else if(input$modelF2=="M2"){
              print(c(input$modelF1,input$modelF2))
              temp.result<-TwoComp_Volume_RateConstant(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],sel.data[i,5])
              print(temp.result)
            } else if(input$modelF2=="M3"){
              temp.result<-TwoComp_Volume_Clearance_HalfLife(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],sel.data[i,5])
            } else if(input$modelF2=="M4"){
              temp.result<-TwoComp_Coefficient_Exponent(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],sel.data[i,5])
            } else if(input$modelF2=="M5"){
              temp.result<-TwoComp_Volume_Exponent(sel.data[i,2],sel.data[i,3],
                                    sel.data[i,4],sel.data[i,5])
            }
            tot.result<-rbind(tot.result,c(i,temp.result[,2]))
            if(i==1) colnames(tot.result)<-c("ID",rownames(temp.result))
          }
        } else if(input$modelF1=="three"){
          for(i in 1:nrow(sel.data)){
            if(input$modelF2=="M1"){
              temp.result<-ThreeComp_Volume_Clearance(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],
                                    sel.data[i,5],sel.data[i,6],sel.data[i,7])
            } else if(input$modelF2=="M2"){
              temp.result<-ThreeComp_Volume_RateConstant(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],
                                    sel.data[i,5],sel.data[i,6],sel.data[i,7])
            } else if(input$modelF2=="M3"){
              temp.result<-ThreeComp_Volume_Clearance_HalfLife(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],
                                    sel.data[i,5],sel.data[i,6],sel.data[i,7])
            } else if(input$modelF2=="M4"){
              temp.result<-ThreeComp_Coefficient_Exponent(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],
                                    sel.data[i,5],sel.data[i,6],sel.data[i,7])
            } else if(input$modelF2=="M5"){
              temp.result<-ThreeComp_Volume_Exponent(sel.data[i,2],
                                    sel.data[i,3],sel.data[i,4],
                                    sel.data[i,5],sel.data[i,6],sel.data[i,7])
            }
            tot.result<-rbind(tot.result,c(i,temp.result[,2]))
            if(i==1) colnames(tot.result)<-c("ID",rownames(temp.result))}}
        out<-tot.result
        PKconverter.env$tot.result<-tot.result})

      output$SaveIndivResult<-shiny::downloadHandler(
        filename = function() {
          paste("IndivParam.csv")},
        content = function(file){
          utils::write.table(PKconverter.env$tot.result, file, sep = ",",
                      row.names = FALSE,col.names = TRUE)})
      }#end server
    )#end App
    shiny::runApp(PK_app,launch.browser=TRUE)
}
