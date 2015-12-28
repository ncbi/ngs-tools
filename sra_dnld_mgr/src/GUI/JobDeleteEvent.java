package GUI;

import javax.swing.JOptionPane;
import job.JobData;

public class JobDeleteEvent
{
    private final MainWindow main_window;
    
    public void delete_job( JobData job )
    {
        int response = JOptionPane.showConfirmDialog(
                main_window,
                String.format( "delete job %s", job.get_short_source() ),
                "question",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.INFORMATION_MESSAGE );
        
        if ( response == JOptionPane.YES_OPTION )
            main_window.RemoveJobAndReArrange( job );
    }
    
    JobDeleteEvent( MainWindow main_window )
    {
        this.main_window = main_window;
    }
}