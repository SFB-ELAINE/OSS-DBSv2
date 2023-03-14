from ossdbs.Output import Output


class TestOutput:

    def test_output_path(self):
        output = Output(path='testresult/er')
        output.output_path()
