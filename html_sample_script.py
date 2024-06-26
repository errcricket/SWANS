import os
import sys

project = sys.argv[1]

from sample_list import get_samples
SAMPLE_LIST = get_samples(project)

root = project + '/'

for s in SAMPLE_LIST:
    filename = 'html_report_' + project + '/' + s + '.Rmd'
    path = '../data/endpoints/' + project  + '/' + s + '/'

    with open(filename, 'w') as output_file:
        newline = '\n'
        header = '---\ntitle: \''+ s +'\'\noutput: html_document\nheader-includes:\n  - \\usepackage{pdfpages}\n--- \n'
        link = '[Web summary](' + path + '10X/outs/web_summary.html)\n'
        setup = '```{r include = FALSE}\nlibrary(magick)\nlibrary(knitr)\n\n'
        sample_name = 'sample_name <- \'' + s + '\'\n'
        create_imdir = 'dir.create(paste0(sample_name, \'_images\'))\n'
        sweep_pkval = 'pk_val<-read.csv(\''+ path + 'tables/' + s + '_' + project + '_pk_value.txt\', header=FALSE)\n'
        sweep_header = '### Sweep Plot & Table \n'
        display_pkval = 'The peak value is `r pk_val[1,1]`\n\n'
        imsetup = '```{r echo = FALSE, out.height="50%", out.width="80%"}\n'
        swplot = 'sweep_plot<-image_read_pdf(\'' + path + 'figures/pk_sweep_plot_' + s + '_' + project + '.pdf\')\n'
        export_swplot = 'sweep_plot[1] %>% image_write(., path = paste0(sample_name, \'_images/\', sample_name, \'_sweepimage.png\'), format = \'png\')\n'
        show_swplot = 'knitr::include_graphics(paste0(sample_name, \'_images/\', sample_name, \'_sweepimage.png\'))\n'
        swtable = 'sweep_table<-read.table(\''+ path + 'tables/' + s + '_' + project + '_pk_values.txt\', header=TRUE, sep=\'\\t\')\nkable(sweep_table)\n' 

        db_header = '### Doublet Plots \n'
        dbplot = 'doublet_pages<-image_read_pdf(\'' + path + 'figures/' + s + '_' + project + '_doublet.pdf\')\n'
        init_j = 'j <- 1:2\n'
        loop_exportdb = 'for(i in j) {\ndoublet_pages[i] %>% image_write(., path = paste0(sample_name, \'_images/\', sample_name, \'_image\',i,\'.png\'), format = \'png\')\n}\n'
        show_db1 = 'knitr::include_graphics(paste0(sample_name, \'_images/\', sample_name, \'_image1.png\'))\n'
        show_db2 = 'knitr::include_graphics(paste0(sample_name, \'_images/\', sample_name, \'_image2.png\'))\n'
        end_codechunk = '```\n'

        output_file.write(header)
        output_file.write(link)
        
        output_file.write(setup)
        output_file.write(sample_name)
        output_file.write(create_imdir)
        output_file.write(sweep_pkval)
        output_file.write(end_codechunk)
        
        output_file.write(sweep_header)
        output_file.write(display_pkval)
        output_file.write(imsetup)
        output_file.write(swplot)
        output_file.write(export_swplot)
        output_file.write(show_swplot)
        output_file.write(swtable)
        output_file.write(end_codechunk)
        
        output_file.write(db_header)
        output_file.write(imsetup)
        output_file.write(dbplot)
        output_file.write(init_j)
        output_file.write(loop_exportdb)
        output_file.write(show_db1)
        output_file.write(show_db2)
        output_file.write(end_codechunk)
