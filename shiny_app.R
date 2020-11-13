library(shiny)
library(plotly)
library(dplyr)
library(reshape2)
library(shinyFiles)


# Define UI for application that draws a histogram
ui = fluidPage(
  
  # Application title
  titlePanel("IDR Visualization"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      textInput(
        inputId = "idrStart", value = '', label = "IDR Start Residue",
        placeholder = "e.g. 1"
      ),
      textInput(
        inputId = "idrEnd", value = '', label = "IDR End Residue",
        placeholder = "e.g. 31"
      ),
      textInput(
        inputId = "fdStart", value = '', label = "FD Start Residue",
        placeholder = "e.g. 32 or enter NA if no FD"
      ),
      textInput(
        inputId = "fdEnd", value = '', label = "FD End Residue",
        placeholder = "e.g. 200 or enter NA if no FD"
      ),
      textInput(
        inputId = "num_clusters", value = '', label = "Num Clusters",
        placeholder = "# of Clusters (3)"
      ),
      textInput(
        inputId = "rmsdresolution", value = '', label = "RMSD Resolution",
        placeholder = "# of CA atoms to calc RMSD (10)"
      ),
      textInput(
        inputId = "framestride", value = '', label = "Frame Stride",
        placeholder = "Use every nth frame (1)"
      ),
      checkboxInput(inputId = 'show_mean_cluster', label = 'Display Mean Trajectory'),
      shinyUI(bootstrapPage(
        shinyFilesButton('pdbfile', 'File to PDB', 'Select File to PDB', multiple=FALSE),
        verbatimTextOutput("pdbfilechosen", placeholder = TRUE)
      )),
      shinyUI(bootstrapPage(
        shinyFilesButton('trajfile', 'File to Traj', 'Select File to Traj File', multiple=FALSE),
        verbatimTextOutput("trajfilechosen", placeholder = TRUE)
      )),
      shinyUI(bootstrapPage(
        shinyDirButton('trajfolder', 'Folder to Traj', 'Select Parent Folder Storing Traj Files'),
        verbatimTextOutput("trajfolderchosen", placeholder = TRUE)
      ))
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      plotlyOutput('plot')
    )
  )
)



# Define server logic required to draw a histogram
shinyServer = function(input, output) {
  
  output$plot <- renderPlotly({
    
    idr_start = as.integer(input$idrStart)
    idr_end = as.integer(input$idrEnd)
    fd_start = as.integer(input$fdStart)
    fd_end = as.integer(input$fdEnd)
    num_clusters = as.integer(input$num_clusters)
    rmsd_resolution = as.integer(input$rmsdresolution)
    frame_stride = as.integer(input$framestride)
    
    fd_start_bool = (!is.na(fd_start) & fd_start != '') | is.na(fd_start) 
    fd_end_bool = (!is.na(fd_end) & fd_end != '') | is.na(fd_end) 
    
    traj_present =  global_traj_file$datapath != getwd() | global_traj_folder$datapath != getwd()
    
    if (idr_start != '' & idr_end != '' & fd_start_bool & fd_end_bool & num_clusters != '' & rmsd_resolution != '' & frame_stride != '' & global_pdb_file$datapath != getwd() & traj_present) {
      #run python script
      print("Generating Visualization Input")
      
      if (global_traj_file$datapath != getwd() & global_traj_folder$datapath != getwd()) {
        traj_input = ''
      } else if (global_traj_file$datapath != getwd()) {
        traj_input = global_traj_file$datapath
      } else {
        traj_input = global_traj_folder$datapath
      }
      
      python_cmd = paste('/Users/ishan/opt/anaconda3/bin/python', './gen_structural_clusters.py', '-s', idr_start, '-e', idr_end, '-t', fd_start, '-f', fd_end, '-p', global_pdb_file$datapath,  '-j', traj_input, '-r', rmsd_resolution, '-n', num_clusters, '-d', frame_stride, sep = ' ')
      system(python_cmd)
      
      #plot output
      print("Plotting Output")
      pdb_info = read.csv('./pdb_info.csv')
      pdb_output = tryCatch(read.csv('./pdb_coordinate_representation.csv'), error = function(e) return(data.frame()))
      cluster_output = read.csv('./cluster_output.csv')
      
      pdb_info$text = paste(pdb_info$residue_num, pdb_info$residue_name, pdb_info$atom, sep = '_')
      
      
      if (input$show_mean_cluster == FALSE)
      {
        cluster_output_wide = data.frame()
        
        cluster_output_wide_x = dcast(cluster_output, aa_num + frame_num ~ cluster, value.var='x')
        cluster_output_wide_y = dcast(cluster_output, aa_num + frame_num ~ cluster, value.var='y')
        cluster_output_wide_z = dcast(cluster_output, aa_num + frame_num ~ cluster, value.var='z')
        cluster_output_wide_frame_num = dcast(cluster_output, aa_num + frame_num ~ cluster, value.var='frame_num')
        
        for (i in 3:ncol(cluster_output_wide_x))
        {
          curr_colname = colnames(cluster_output_wide_x)[i]
          colnames(cluster_output_wide_x)[i] = paste('x', curr_colname, sep = '_')
          colnames(cluster_output_wide_y)[i] = paste('y', curr_colname, sep = '_')
          colnames(cluster_output_wide_z)[i] = paste('z', curr_colname, sep = '_')
          colnames(cluster_output_wide_frame_num)[i] = paste('frame_num', curr_colname, sep = '_')
        }
        
        output_wide = left_join(cluster_output_wide_x, cluster_output_wide_y, by = c('aa_num', 'frame_num')) %>% left_join(cluster_output_wide_z, by = c('aa_num', 'frame_num')) %>% left_join(cluster_output_wide_frame_num, by = c('aa_num', 'frame_num')) %>% as.data.frame()
        
        if (nrow(pdb_output) > 0)
        {
          output_wide$pdb_x = NA
          output_wide$pdb_y = NA
          output_wide$pdb_z = NA
          output_wide$text = NA
          output_wide$charge = NA
          
          num_na_rows_to_add = nrow(pdb_output) - nrow(output_wide)
          if (num_na_rows_to_add > 0)
          {
            output_wide[(nrow(output_wide)+1):(nrow(output_wide)+num_na_rows_to_add),] = NA
          }
          
          output_wide$pdb_x[1:nrow(pdb_output)] = pdb_output$x
          output_wide$pdb_y[1:nrow(pdb_output)] = pdb_output$y
          output_wide$pdb_z[1:nrow(pdb_output)] = pdb_output$z
          output_wide$text[1:nrow(pdb_info)] = pdb_info$text
          output_wide$charge[1:nrow(pdb_info)] = pdb_info$charge
          
          output_wide$charge = as.factor(output_wide$charge)
        }
        
        p = plot_ly()
        unique_clusters = unique(output_wide$cluster)
        color_list = c('#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52')
        for (i in 1:length(unique(cluster_output$cluster)))
        {
          curr_cluster = unique(cluster_output$cluster)[i]
          p = add_trace(p, x = output_wide[[paste('x', curr_cluster, sep = '_')]], y = output_wide[[paste('y', curr_cluster, sep = '_')]], z = output_wide[[paste('z', curr_cluster, sep = '_')]], type = 'scatter3d', mode = 'lines+markers', color = factor(output_wide[[paste('frame_num', curr_cluster, sep = '_')]]),
                        opacity = 1, line = list(color=color_list[i], width = 1),
                        legendgroup = curr_cluster, 
                        marker = list(color=color_list[i], size = 2))
        }
        if (nrow(pdb_output) > 0)
        {
          p = add_trace(p, x = output_wide[['pdb_x']], y = output_wide[['pdb_y']], z = output_wide[['pdb_z']], type = 'scatter3d', mode = 'markers', 
                        marker = list(color = factor(output_wide[['charge']],labels=c("red","grey","blue")), size = 4),
                        text = output_wide[['text']], hoverinfo = 'text')
        }
        p = p %>% layout(showlegend = TRUE)
        p 
      }
      else
      {
        cluster_output_avg = cluster_output %>% group_by(aa_num, cluster) %>% summarise(mean_x = mean(x), sd_x = sd(x), mean_y = mean(y), sd_y = sd(y), mean_z = mean(z), sd_z = sd(z), n = n()) %>% as.data.frame()
        cluster_output_avg = cluster_output_avg %>% mutate(se_x = sd_x/sqrt(n), se_y = sd_y/sqrt(n), se_z = sd_z/sqrt(n)) %>% as.data.frame()
        cluster_output_avg = cluster_output_avg %>% mutate(sd_x_upper = mean_x+1.96*sd_x, 
                                                           sd_x_lower = mean_x-1.96*sd_x,
                                                           sd_y_upper = mean_y+1.96*sd_y, 
                                                           sd_y_lower = mean_y-1.96*sd_y,
                                                           sd_z_upper = mean_z+1.96*sd_z, 
                                                           sd_z_lower = mean_z-1.96*sd_z) %>% as.data.frame()
        
        
        cluster_output_wide_avg = data.frame()
        
        cluster_output_wide_avg_x_m = dcast(cluster_output_avg, aa_num ~ cluster, value.var='mean_x')
        cluster_output_wide_avg_y_m = dcast(cluster_output_avg, aa_num ~ cluster, value.var='mean_y')
        cluster_output_wide_avg_z_m = dcast(cluster_output_avg, aa_num ~ cluster, value.var='mean_z')
        
        cluster_output_wide_avg_x_sd_u = dcast(cluster_output_avg, aa_num ~ cluster, value.var='sd_x_upper')
        cluster_output_wide_avg_x_sd_l = dcast(cluster_output_avg, aa_num ~ cluster, value.var='sd_x_lower')
        cluster_output_wide_avg_y_sd_u = dcast(cluster_output_avg, aa_num ~ cluster, value.var='sd_y_upper')
        cluster_output_wide_avg_y_sd_l = dcast(cluster_output_avg, aa_num ~ cluster, value.var='sd_y_lower')
        cluster_output_wide_avg_z_sd_u = dcast(cluster_output_avg, aa_num ~ cluster, value.var='sd_z_upper')
        cluster_output_wide_avg_z_sd_l = dcast(cluster_output_avg, aa_num ~ cluster, value.var='sd_z_lower')
        
        for (i in 2:ncol(cluster_output_wide_avg_x_m))
        {
          curr_colname = colnames(cluster_output_wide_avg_x_m)[i]
          colnames(cluster_output_wide_avg_x_m)[i] = paste('mean_x', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_y_m)[i] = paste('mean_y', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_z_m)[i] = paste('mean_z', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_x_sd_u)[i] = paste('sd_x_upper', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_x_sd_l)[i] = paste('sd_x_lower', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_y_sd_u)[i] = paste('sd_y_upper', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_y_sd_l)[i] = paste('sd_y_lower', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_z_sd_u)[i] = paste('sd_z_upper', curr_colname, sep = '_')
          colnames(cluster_output_wide_avg_z_sd_l)[i] = paste('sd_z_lower', curr_colname, sep = '_')
        }
        
        output_avg_wide = left_join(cluster_output_wide_avg_x_m, cluster_output_wide_avg_y_m, by = c('aa_num')) %>% 
          left_join(cluster_output_wide_avg_z_m, by = c('aa_num')) %>% left_join(cluster_output_wide_avg_x_sd_u, by = c('aa_num')) %>% 
          left_join(cluster_output_wide_avg_x_sd_l, by = c('aa_num')) %>% left_join(cluster_output_wide_avg_y_sd_u, by = c('aa_num')) %>% 
          left_join(cluster_output_wide_avg_y_sd_l, by = c('aa_num')) %>% left_join(cluster_output_wide_avg_z_sd_u, by = c('aa_num')) %>% 
          left_join(cluster_output_wide_avg_z_sd_l, by = c('aa_num')) %>% as.data.frame()
        
        if (nrow(pdb_output) > 0)
        {
          output_avg_wide$pdb_x = NA
          output_avg_wide$pdb_y = NA
          output_avg_wide$pdb_z = NA
          output_avg_wide$text = NA
          output_avg_wide$charge = NA
          
          num_na_rows_to_add = nrow(pdb_output) - nrow(output_avg_wide)
          if (num_na_rows_to_add > 0)
          {
            output_avg_wide[(nrow(output_avg_wide)+1):(nrow(output_avg_wide)+num_na_rows_to_add),] = NA
          }
          
          output_avg_wide$pdb_x[1:nrow(pdb_output)] = pdb_output$x
          output_avg_wide$pdb_y[1:nrow(pdb_output)] = pdb_output$y
          output_avg_wide$pdb_z[1:nrow(pdb_output)] = pdb_output$z
          output_avg_wide$text[1:nrow(pdb_info)] = pdb_info$text
          output_avg_wide$charge[1:nrow(pdb_info)] = pdb_info$charge
          
          output_avg_wide$charge = as.factor(output_avg_wide$charge)
        }
        
        p = plot_ly()
        unique_clusters = unique(cluster_output_avg$cluster)
        color_list = c('#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52')
        
        for (i in 1:length(unique(cluster_output_avg$cluster)))
        {
          curr_cluster = unique(cluster_output_avg$cluster)[i]
          n = cluster_output_avg[cluster_output_avg$cluster == curr_cluster,'n'][1]
          
          if (n >= 25)
          {
            p = add_trace(p, x = output_avg_wide[[paste('mean_x', curr_cluster, sep = '_')]], y = output_avg_wide[[paste('mean_y', curr_cluster, sep = '_')]], z = output_avg_wide[[paste('mean_z', curr_cluster, sep = '_')]], type = 'scatter3d', mode = 'lines+markers',
                          opacity = 1, line = list(color=color_list[i], width = 4),
                          legendgroup = curr_cluster,
                          marker = list(color=color_list[i], size = 3))
            p = add_trace(p, x = output_avg_wide[[paste('sd_x_lower', curr_cluster, sep = '_')]], y = output_avg_wide[[paste('sd_y_lower', curr_cluster, sep = '_')]], z = output_avg_wide[[paste('sd_z_lower', curr_cluster, sep = '_')]], type = 'scatter3d', mode = 'lines',
                          opacity = .5, line = list(color=color_list[i], width = 2),
                          legendgroup = curr_cluster)
            p = add_trace(p, x = output_avg_wide[[paste('sd_x_upper', curr_cluster, sep = '_')]], y = output_avg_wide[[paste('sd_y_upper', curr_cluster, sep = '_')]], z = output_avg_wide[[paste('sd_z_upper', curr_cluster, sep = '_')]], type = 'scatter3d', mode = 'lines',
                          opacity = .5, line = list(color=color_list[i], width = 2),
                          legendgroup = curr_cluster)
          }
        }
        if (nrow(pdb_output) > 0)
        {
          p = add_trace(p, x = output_avg_wide[['pdb_x']], y = output_avg_wide[['pdb_y']], z = output_avg_wide[['pdb_z']], type = 'scatter3d', mode = 'markers', 
                        marker = list(color = factor(output_avg_wide[['charge']],labels=c("red","grey","blue")), size = 4),
                        text = output_avg_wide[['text']], hoverinfo = 'text')
        }
        p = p %>% layout(showlegend = TRUE)
        p
        
      }
    }
    
  })
  
  
  shinyFileChoose(input, 'pdbfile', roots = c(home = '~'), filetypes=c('', 'pdb'))
  shinyFileChoose(input, 'trajfile', roots = c(home = '~'), filetypes=c('', 'xtc', 'dcd'))
  shinyDirChoose(input, 'trajfolder', roots = c(home = '~'), filetypes=c(''))
  
   pdbf = reactive(input$pdbfile)
   output$pdbfilechosen <- renderText({
     as.character(parseFilePaths(c(home = "~"),pdbf())$datapath)
   })
   
   global_pdb_file = reactiveValues(datapath = getwd())
   observeEvent(input$pdbfile, {
     global_pdb_file$datapath = as.character(parseFilePaths(c(home = "~"),pdbf())$datapath)
   })
   
   trajf = reactive(input$trajfile)
   output$trajfilechosen <- renderText({
     as.character(parseFilePaths(c(home = "~"),trajf())$datapath)
   })
   
   global_traj_file = reactiveValues(datapath = getwd())
   observeEvent(input$trajfile, {
     global_traj_file$datapath = as.character(parseFilePaths(c(home = "~"),trajf())$datapath)
   })
   

  global_traj_folder = reactiveValues(datapath = getwd())

  trajfolder = reactive(input$trajfolder)
  output$trajfolder = renderText({
     global_traj_folder$datapath
  })
  observeEvent(ignoreNULL = TRUE,
              eventExpr = {
                input$trajfolder
              },
              handlerExpr = {
                if (!"path" %in% names(trajfolder())) return()
                home = normalizePath("~")
                global_traj_folder$datapath = file.path(home, paste(unlist(trajfolder()$path[-1]), collapse = .Platform$file.sep))
              })
  
}


# Run the application
shinyApp(ui = ui, server = shinyServer)