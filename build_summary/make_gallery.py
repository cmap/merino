import jinja2


def mk_gal(images, outfolder):

    templateLoader = jinja2.FileSystemLoader(searchpath="/Users/elemire/Workspace/merino/build_summary/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "gallery_lazy.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    opt = {'ncol':4, 'width': 230, 'height': 225, 'summary': ''}
    outputText = template.render(imlist=images, opt=opt)  # this is where to put args to the template renderer




    #print(outputText)

    with open(outfolder, "wb") as fh:
        fh.write(outputText)