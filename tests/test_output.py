from ossdbs.output import OutputDirectory


class TestOutput:

    def test_output_path(self):
        output = OutputDirectory(directory='testresult/subfolder/subfolder')
        # output.directory()
