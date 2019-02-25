import os
import jinja2

path = os.getcwd()
templateLoader = jinja2.FileSystemLoader(searchpath=path)
templateEnv = jinja2.Environment(loader=templateLoader)

def mk_gal(images, outfolder):

    TEMPLATE_FILE = "gallery_lazy.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    opt = {'ncol':4, 'width': 230, 'height': 225, 'summary': ''}
    outputText = template.render(imlist=images, opt=opt)  # this is where to put args to the template renderer

    #print(outputText)

    with open(outfolder, "wb") as fh:
        fh.write(outputText)

def mk_index(table_headers, table_tuples, outfolder, project_name):

    TEMPLATE_FILE = "table_template.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    opt = {'title':project_name}
    outputText = template.render(headers=table_headers, summary=table_tuples, opt=opt)

    with open(os.path.join(outfolder, "full_gallery.html"), "wb") as fh:
        fh.write(outputText)