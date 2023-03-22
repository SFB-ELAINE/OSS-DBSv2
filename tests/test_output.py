from ossdbs.outputX import OutputDirectory


class TestOutput:

    def test_output_path(self):
        output = OutputDirectory(directory='testresult/er')
        output.output_directory()
