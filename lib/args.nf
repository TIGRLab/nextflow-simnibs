nextflow.enable.dsl = 2

def validateArgs(required){
/*
* Check if required arguments are missing
*/
     required.grep{(it.value == null || it.value == "")}
}

def getUsage(usage_file, bindings){
/*
* Make usage message
*/
    def usage = new File(usage_file)

    groovy.text.SimpleTemplateEngine()
          .createTemplate(usage.text)
          .make(bindings)
}
