import groovy.text.SimpleTemplateEngine;

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
    def engine = new SimpleTemplateEngine();
    engine.createTemplate(usage.text).make(bindings)
}
