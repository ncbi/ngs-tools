package data;

import java.util.concurrent.ConcurrentLinkedQueue;

public class QConnect< EntryType >
{
    private final ConcurrentLinkedQueue< EntryType > active;
    private final ConcurrentLinkedQueue< EntryType > stock;
    private final int max_created;
    private int created;
    
    public EntryType make_entry() { return null; }
    
    public synchronized EntryType get_from_stock()
    {
        EntryType res = stock.poll();
        if ( res == null && created < max_created )
        {
            res = make_entry();
            if ( res != null ) created++;
        }
        return res;
    }

    public synchronized void put_to_stock( final EntryType e ) {  stock.offer( e ); }
    public synchronized boolean has_entries() { return !active.isEmpty(); }
    public synchronized void put( final EntryType e ) { active.offer( e ); }
    public synchronized EntryType get() { return active.poll(); }
    
    public synchronized void clear()
    {
        active.clear();
        stock.clear();
        created = 0;
    }
    
    public QConnect( final int max_created )
    {
        active = new ConcurrentLinkedQueue<EntryType>();
        stock = new ConcurrentLinkedQueue<EntryType>();
        created = 0;
        this.max_created = max_created;
    }
    
}
