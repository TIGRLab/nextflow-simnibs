import groovy.text.SimpleTemplateEngine;
import java.util.LinkedHashMap;
import org.apache.commons.lang.WordUtils;

def validateArgs(LinkedHashMap required){
/*
* Check if required arguments are missing
*/

     required.grep { (it.value == null || it.value == "") }
}


@groovy.transform.TupleConstructor class Argument {
    /*
    * Simple data class for representing Nextflow command-line arguments
    */

    String flag, description, metavar, value
    Boolean optional

    public String makeEntry() {
        /*
        * Construct argument usage
        */
        "${this.flag}\t${this.description}".toString()
    }

    public String makeUsage(){
        "${this.flag} ${(this.metavar ? '<' + this.metavar + '>' : '')}".toString()
    }

    public Boolean validArg(){
        if (!this.optional){
            return !this.value.isEmpty()
        }

        true
    }

 }

@groovy.transform.TupleConstructor class ArgumentParser {
    /*
    * Rough implementation of Python-style API for argument parsing
    */

    List required = [], optional = []
    String title, description, note, scriptName

    ArgumentParser(String title,
    String description,
    String scriptName,
    String note=""){
        this.title = title
        this.description = description
        this.scriptName = scriptName
        this.note = note
    }
    
    public void addArgument(String flag, String description, String value, String metavar=null) {
        this.required.push(
            new Argument(
            flag: flag,
            description: description,
            metavar: metavar,
            value: value.toString(),
            optional: false
            )
        )
    }

    public void addOptional(String flag, String description, String metavar=null) {
        this.optional.push(
            new Argument(
            flag: flag,
            description: description,
            metavar: metavar,
            value: null,
            optional: true
            )
        )
    }

    private String generateArgs() {
        /*
        * Generate required section
        */

        def formatted = this.required
            .collect { arg -> "|\t" + arg.makeEntry() }
            .join("\n")
        """\

        |REQUIRED
        ${formatted}
        |\t-c\tPath to nextflow config file
        """.toString().stripIndent().stripMargin()
    }

    private String generateOpts() {
        /*
        * Generate optional section
        */

        def formatted = this.optional
            .collect { arg -> "|\t" + arg.makeEntry() }
            .join("\n")
        """\

        |OPTIONAL
        ${formatted}
        |\t--help\tPrint this usage log
        """.toString().stripIndent().stripMargin()
    }


    private String generateUsage() {
        /*
        * Generate usage doc
        */

        def usage_args = this.required
            .collect { arg -> arg.makeUsage() }
            .join(" ")

        def usage = """\
        |USAGE

        |\tnextflow run ${this.scriptName} ${(!this.optional.isEmpty() ? '[options...]' : '')} -c <CONFIG> ${usage_args}
        """.toString().stripIndent().stripMargin()

        WordUtils.wrap(usage, 100)
    }

    public void displayDoc() {
        def headerStr = """\
        ${this.title}
        ${'=' * this.title.size()}

        ${description}
        """.stripIndent()
        println(headerStr)
        println(this.generateUsage())
        println(this.generateArgs())
        println(this.generateOpts())
        println(this.note.stripIndent())
    }

    public String makeDoc() {
        /*
        * Make usage documentation
        */


        """\
        |${this.title}
        |${'=' * this.title.size()}

        |${description}

        |${this.generateUsage()}
        |${this.generateArgs()}
        |${this.generateOpts()}

        |${(this.note ? 'NOTE' : '')}
        |${this.note}
        """.stripIndent().stripMargin()
    }

    public List isMissingRequired(){
        /*
        * Return all args with invalid values
        */
        this.required
            .findAll{ arg -> !arg.validArg() }
            .collect{ arg -> arg.flag }
    }

 }


def getArgumentParser(Map args){
    /*
    * This is necessary for Nextflow 'include' to work
    * Classes cannot be directly imported so we have
    * to wrap their construction with a very simple
    * constructor function
    *
    * Arguments:
    *   args (Map): 
    *       title (String): Title of program
    *       description (String): Description text
    *       note (String): Additional notes
    *
    * Returns:
    *   argumentParser (ArgumentParser): Class for handling arguments
    */

    def arg = new ArgumentParser(args.get('title'),
        args.get('description'),
        args.get('scriptName'),
        args.get('note')
    )
    arg
}
