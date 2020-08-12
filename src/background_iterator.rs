use std::sync::mpsc::{Receiver, sync_channel};
use std;
use std::thread;

/// Execute an iterator on a worker thread, which can work ahead a configurable number of items.
pub struct BackgroundIterator<T> {
    rx: Receiver<Option<T>>,
    done: bool,
}

impl<T: Send> Iterator for BackgroundIterator<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        if self.done {
            return None;
        }

        match self.rx.recv() {
            Ok(Some(v)) => Some(v),
            Ok(None) => {
                self.done = true;
                None
            }

            // if the producer thread dies, the iterator
            // silently exits, without an error.
            Err(_) => None,
        }
    }
}

impl<T: 'static + Send> BackgroundIterator<T> {
    /// Iterate through `itr` on a newly created thread, and send items back to the returned
    /// `BackgroundIterator` for consumption on the calling thread. The worker thread will
    /// continue to produce elements until it is `max_read_ahead` items ahead of the consumer iterator.
    pub fn new<I: 'static + Send + Iterator<Item = T>>(
        itr: I,
        max_read_ahead: usize,
    ) -> BackgroundIterator<T> {
        let (tx, rx) = sync_channel::<Option<T>>(max_read_ahead);
        let _ = thread::spawn(move || {
            for item in itr {
                if tx.send(Some(item)).is_err() {
                    return;
                };
            }

            tx.send(None).unwrap();
        });

        BackgroundIterator { rx, done: false }
    }
}
